const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});

    // Domyślnie ReleaseFast, bo to projekt pomiarowy: benchmark w trybie
    // Debug mierzyłby instrumentację bezpieczeństwa, a nie algorytm.
    // Można nadpisać: zig build bench -Doptimize=Debug
    const optimize = b.option(
        std.builtin.OptimizeMode,
        "optimize",
        "Tryb optymalizacji (domyślnie ReleaseFast)",
    ) orelse .ReleaseFast;

    const mrs_mod = b.addModule("mrs", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const exe = b.addExecutable(.{
        .name = "mrs-lab",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "mrs", .module = mrs_mod },
            },
        }),
    });
    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| run_cmd.addArgs(args);
    const run_step = b.step("run", "Uruchom mrs (domyślnie: demo)");
    run_step.dependOn(&run_cmd.step);

    // --- testy -------------------------------------------------------------
    // Testy uruchamiamy w TRZECH trybach, nie w jednym. Powód: domyślny tryb
    // tego projektu to ReleaseFast (bo to projekt pomiarowy), a w ReleaseFast
    // kontrole bezpieczeństwa są WYŁĄCZONE i przepełnienie czy wyjście poza
    // zakres są cichym UB. Tryb Debug i ReleaseSafe zamieniają to na panikę,
    // więc dopiero one sprawdzają stabilność. ReleaseFast zostaje, żeby
    // złapać błędy zależne od optymalizacji.
    const test_step = b.step("test", "Uruchom wszystkie testy (Debug + ReleaseSafe + ReleaseFast)");

    const test_modes = [_]std.builtin.OptimizeMode{ .Debug, .ReleaseSafe, .ReleaseFast };
    for (test_modes) |mode| {
        const lib_mod = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
            .optimize = mode,
        });
        const lib_tests = b.addTest(.{ .root_module = lib_mod });
        test_step.dependOn(&b.addRunArtifact(lib_tests).step);

        const exe_mod = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = mode,
            .imports = &.{
                .{ .name = "mrs", .module = lib_mod },
            },
        });
        const exe_tests = b.addTest(.{ .root_module = exe_mod });
        test_step.dependOn(&b.addRunArtifact(exe_tests).step);
    }

    // --- komendy wygodne ---------------------------------------------------
    const demo_step = b.step("demo", "Pokaz matematyczny MRS-0");
    const demo_cmd = b.addRunArtifact(exe);
    demo_cmd.step.dependOn(b.getInstallStep());
    demo_cmd.addArg("demo");
    demo_step.dependOn(&demo_cmd.step);

    const verify_step = b.step("verify", "Sprawdzenia spójności (kod wyjścia 1 przy błędzie)");
    const verify_cmd = b.addRunArtifact(exe);
    verify_cmd.step.dependOn(b.getInstallStep());
    verify_cmd.addArg("verify");
    verify_step.dependOn(&verify_cmd.step);

    const bench_step = b.step("bench", "Raport pomiarowy T1–T4 do results/RESULTS.md");
    const bench_cmd = b.addRunArtifact(exe);
    bench_cmd.step.dependOn(b.getInstallStep());
    bench_cmd.addArgs(&.{ "bench", "--out", "results/RESULTS.md" });
    bench_step.dependOn(&bench_cmd.step);

    // Benchmark smoke run in a CHECKED mode. The benchmark binary itself is
    // built ReleaseFast (that is the point of a benchmark), which means a bounds
    // error inside benchmark code is a segfault rather than a panic. This step
    // builds the same code ReleaseSafe and runs a short series, so that class of
    // defect surfaces as a readable panic. It exists because a real bounds bug in
    // a precomputation path was found exactly this way.
    const check_step = b.step("bench-check", "Short benchmark run built ReleaseSafe (catches bounds errors)");
    const check_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = .ReleaseSafe,
        .imports = &.{
            .{ .name = "mrs", .module = b.createModule(.{
                .root_source_file = b.path("src/root.zig"),
                .target = target,
                .optimize = .ReleaseSafe,
            }) },
        },
    });
    const check_exe = b.addExecutable(.{ .name = "mrs-lab-check", .root_module = check_mod });
    const check_cmd = b.addRunArtifact(check_exe);
    check_cmd.addArgs(&.{ "bench", "--quick" });
    check_step.dependOn(&check_cmd.step);

    const explore_step = b.step("explore", "Tabela wyników P2/P4 (do sumy p+q+r <= 4)");
    const explore_cmd = b.addRunArtifact(exe);
    explore_cmd.step.dependOn(b.getInstallStep());
    explore_cmd.addArg("explore");
    explore_step.dependOn(&explore_cmd.step);

    const explore_full_step = b.step("explore-full", "Tabela wyników P2/P4 dla sumy <= 5 (wolniejsze)");
    const explore_full_cmd = b.addRunArtifact(exe);
    explore_full_cmd.step.dependOn(b.getInstallStep());
    explore_full_cmd.addArgs(&.{ "explore", "--max", "5" });
    explore_full_step.dependOn(&explore_full_cmd.step);

    const bench_quick_step = b.step("bench-quick", "Krótki raport pomiarowy");
    const bench_quick_cmd = b.addRunArtifact(exe);
    bench_quick_cmd.step.dependOn(b.getInstallStep());
    bench_quick_cmd.addArgs(&.{ "bench", "--quick" });
    bench_quick_step.dependOn(&bench_quick_cmd.step);
}
