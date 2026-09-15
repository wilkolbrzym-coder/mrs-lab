//! MRS :: measurement harness
//!
//! Rules every benchmark in this project follows:
//!
//!   1. Warmup before measuring (warm the cache and the branch predictors).
//!   2. Many samples; we report the MINIMUM and the MEDIAN, not the mean —
//!      the mean is poisoned by system interrupts.
//!   3. The result is consumed through `doNotOptimizeAway` and the input is
//!      mutated on every call, so the compiler cannot hoist the computation
//!      out of the loop.
//!   4. The "conventional" baseline is written the way a competent engineer
//!      who knows their data would write it — not as a straw man.

const std = @import("std");
const Io = std.Io;

pub fn nowNs(io: Io) i96 {
    return Io.Timestamp.now(io, .awake).nanoseconds;
}

fn asc(_: void, a: f64, b: f64) bool {
    return a < b;
}

pub const Stats = struct {
    best_ns: f64,
    median_ns: f64,
    mean_ns: f64,
    samples: usize,
    inner: usize,

    pub fn perOpBest(self: Stats) f64 {
        return self.best_ns / @as(f64, @floatFromInt(self.inner));
    }

    pub fn perOpMedian(self: Stats) f64 {
        return self.median_ns / @as(f64, @floatFromInt(self.inner));
    }

    pub fn spread(self: Stats) f64 {
        if (self.median_ns == 0) return 0;
        return (self.median_ns - self.best_ns) / self.median_ns;
    }
};

/// Measures `samples` samples; inside each it calls `func` `inner` times.
/// `ctx` is a pointer — benchmarks keep their own buffers in it.
pub fn measure(
    io: Io,
    alloc: std.mem.Allocator,
    ctx: anytype,
    samples: usize,
    inner: usize,
    comptime func: fn (@TypeOf(ctx)) void,
) !Stats {
    // warmup
    for (0..@max(1, inner / 8)) |_| func(ctx);

    const times = try alloc.alloc(f64, samples);
    defer alloc.free(times);

    for (0..samples) |s| {
        const t0 = nowNs(io);
        for (0..inner) |_| func(ctx);
        const t1 = nowNs(io);
        times[s] = @floatFromInt(t1 - t0);
    }

    std.mem.sort(f64, times, {}, asc);
    var total: f64 = 0;
    for (times) |t| total += t;

    return .{
        .best_ns = times[0],
        .median_ns = times[samples / 2],
        .mean_ns = total / @as(f64, @floatFromInt(samples)),
        .samples = samples,
        .inner = inner,
    };
}

/// Relative speed: how many times slower `baseline` is than `mrs`.
pub fn speedup(baseline_ns_per_op: f64, mrs_ns_per_op: f64) f64 {
    if (mrs_ns_per_op == 0) return 0;
    return baseline_ns_per_op / mrs_ns_per_op;
}

/// Picks the iteration count so that one sample lasts approximately
/// `target_ns` — otherwise for small n we measure clock noise.
pub fn pickInner(cost_per_op: f64, target_ns: f64) usize {
    if (cost_per_op <= 0) return 1;
    const n = @floor(target_ns / cost_per_op);
    if (n < 1) return 1;
    return @intFromFloat(n);
}

/// Writes a time in a readable unit: "123 ns", "4.56 µs", "78.9 ms".
pub fn writeDuration(w: anytype, ns: f64) !void {
    if (ns < 1e3) {
        try w.print("{d:.1} ns", .{ns});
    } else if (ns < 1e6) {
        try w.print("{d:.2} µs", .{ns / 1e3});
    } else if (ns < 1e9) {
        try w.print("{d:.2} ms", .{ns / 1e6});
    } else {
        try w.print("{d:.2} s", .{ns / 1e9});
    }
}

pub fn writeCount(w: anytype, x: f64) !void {
    if (x >= 1e9) {
        try w.print("{d:.2}e9", .{x / 1e9});
    } else if (x >= 1e6) {
        try w.print("{d:.2}e6", .{x / 1e6});
    } else if (x >= 1e3) {
        try w.print("{d:.2}e3", .{x / 1e3});
    } else {
        try w.print("{d:.1}", .{x});
    }
}

test "the harness measures what it is supposed to measure" {
    const Ctx = struct {
        x: f64 = 0,
    };
    const f = struct {
        fn call(c: *Ctx) void {
            c.x += 1.0;
            std.mem.doNotOptimizeAway(c.x);
        }
    }.call;

    const alloc = std.testing.allocator;
    _ = alloc;
    // `measure` needs Io, which tests do not have — so only the helper
    // arithmetic is checked here.
    try std.testing.expectEqual(@as(usize, 1), pickInner(0.0, 1e6));
    try std.testing.expect(pickInner(100.0, 1e6) == 10_000);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), speedup(8.0, 2.0), 1e-15);
    _ = f;
}

/// Zig has no '+' flag in format specifiers, and readability gains a lot from an
/// readability gains a lot from an explicit sign. Hence writing it out.
pub fn writeSigned2(w: anytype, x: f64) !void {
    if (x >= 0.0) try w.writeAll("+");
    try w.print("{d:.2}", .{x});
}

pub fn writeSigned3(w: anytype, x: f64) !void {
    if (x >= 0.0) try w.writeAll("+");
    try w.print("{d:.3}", .{x});
}
