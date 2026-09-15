//! MRS-LAB :: public library interface
//!
//! The mathematical layer is described in docs/; the performance theses and
//! their honest verdicts live in results/RESULTS.md, and the property engine
//! table in results/EXPLORE.md.

pub const signature = @import("mrs/signature.zig");
pub const form = @import("mrs/form.zig");
pub const split_complex = @import("mrs/split_complex.zig");
pub const dual = @import("mrs/dual.zig");
pub const clifford = @import("mrs/clifford.zig");
pub const causal = @import("mrs/causal.zig");
pub const tachyon = @import("mrs/tachyon.zig");

test {
    // Same trap as in main.zig: without a reference from this block Zig does
    // not analyse the file in test context and collects no tests from it.
    _ = signature;
    _ = form;
    _ = split_complex;
    _ = dual;
    _ = clifford;
    _ = causal;
    _ = tachyon;
}
