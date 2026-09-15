//! MRS-0 :: Moduł przyczynowości
//!
//! Dwie relacje, które trzeba rozróżniać — i których nierozróżnianie było
//! najpoważniejszym błędem w specyfikacji projektu:
//!
//!   * SEPARACJA (niezorientowana).  u ⊑ v  ⟺  g(v−u, v−u) ≤ 0.
//!     "Różnica jest czasowa albo zerowa". Symetryczna na separacjach
//!     zerowych, więc NIE jest antysymetryczna, więc nie jest porządkiem.
//!
//!   * PRZYCZYNOWOŚĆ (zorientowana).  u ⪯ v  ⟺  (v−u) leży w domkniętym
//!     stożku przyszłości, czyli klasa różnicy to `temporal` albo
//!     `null_like` ORAZ składowa wzdłuż wybranej strzałki czasu ≥ 0.
//!
//! Orientacja wymaga dodatkowej struktury: wyboru strzałki czasu. Przy
//! p = 1 forma wyznacza ją z dokładnością do znaku, przy p ≥ 2 — nie
//! wyznacza jej wcale. To nie jest szczegół implementacyjny, to jest treść
//! Twierdzenia 5.2.
//!
//! ---------------------------------------------------------------------------
//! WYNIK GŁÓWNY MRS-0
//! ---------------------------------------------------------------------------
//! Twierdzenie 5.1 (przechodniość wymaga jednego czasu).
//!   Dla r = 0 i q ≥ 1: ⪯ jest przechodnia  ⟺  p = 1.
//!   (⟸) Stożek przyszłości jest wypukły, a wypukłość stożka daje
//!       przechodniość: a, b w stożku ⟹ a + b w stożku.
//!   (⟹) Dla p ≥ 2 istnieje świadek — patrz `canonicalWitness`.
//!
//! Twierdzenie 5.2 (⪯ jest porządkiem częściowym dokładnie dla sygnatur
//! lorentzowskich).
//!   Zwrotność: zawsze.
//!   Przechodniość: ⟺ p = 1 (przy r = 0).
//!   Antysymetria: ⟺ r = 0 i p ≤ 1.
//!       Dowód (⟸) dla p = 1, r = 0: jeżeli d i −d są w stożku, to
//!       składowa strzałki daje d_arrow = 0, a wtedy g(d) ≥ 0 wymusza
//!       Zerowanie wszystkich składowych przestrzennych, więc d = 0.
//!       Dowód (⟹) dla p ≥ 2: wektor o zerowej składowej strzałki,
//!       niezerowej drugiej składowej czasowej i dostatecznie małej
//!       przestrzennej należy do stożka razem z przeciwieństwem.
//!       Dla r > 0: każdy wektor jądra jest takim świadkiem.
//!   Wniosek: ⪯ jest porządkiem częściowym ⟺ p = 1 ∧ r = 0 ⟺ `isLorentzian()`.
//!   To znaczy, że predykat `isLorentzian` w module sygnatur nie jest
//!   ozdobnikiem — jest warunkiem koniecznym i wystarczającym na to, żeby
//!   przyczynowość była porządkiem.
//!
//! Twierdzenie 5.3 (r > 0 niszczy przyczynowość).
//!   Jeżeli r > 0, to istnieje wektor jądra d ≠ 0 z g(d,d) = 0 i zerową
//!   składową strzałki czasu, czyli u ⪯ u+d ⪯ u przy d ≠ 0. Relacja
//!   przestaje odróżniać punkty. Przyczynowość wymaga r = 0.

const std = @import("std");
const sig = @import("signature.zig");
const form = @import("form.zig");

const Signature = sig.Signature;
const DiagonalForm = form.DiagonalForm;

pub const MAX_DIM: usize = 8;

pub const OrderError = error{
    TimeArrowOutOfRange,
    TimeArrowNotTemporal,
    DimensionTooLarge,
};

/// Świadek naruszenia własności. Wektory trzymamy na stosie, bo moduł
/// przyczynowości bywa używany w gorących pętlach symulacyjnych.
pub const Witness = struct {
    u: [MAX_DIM]f64 = [_]f64{0} ** MAX_DIM,
    w: [MAX_DIM]f64 = [_]f64{0} ** MAX_DIM,
    v: [MAX_DIM]f64 = [_]f64{0} ** MAX_DIM,
    len: usize = 0,

    pub fn uSlice(self: *const Witness) []const f64 {
        return self.u[0..self.len];
    }
    pub fn wSlice(self: *const Witness) []const f64 {
        return self.w[0..self.len];
    }
    pub fn vSlice(self: *const Witness) []const f64 {
        return self.v[0..self.len];
    }

    pub fn writeTo(self: *const Witness, w: anytype) !void {
        try w.writeAll("u=[");
        try writeVec(w, self.uSlice());
        try w.writeAll("]  w=[");
        try writeVec(w, self.wSlice());
        try w.writeAll("]  v=[");
        try writeVec(w, self.vSlice());
        try w.writeAll("]");
    }
};

fn writeVec(w: anytype, v: []const f64) !void {
    for (v, 0..) |x, i| {
        if (i != 0) try w.writeAll(", ");
        try w.print("{d:.3}", .{x});
    }
}

pub const ConeMode = enum {
    /// Klasa `temporal`, skierowany w przyszłość.
    timelike_future,
    /// Klasa `null_like`, skierowany w przyszłość.
    null_future,
};

/// Co się psuje w danej sygnaturze.
pub const Causality = struct {
    reflexive: bool,
    transitive: bool,
    antisymmetric: bool,

    pub fn isPartialOrder(self: Causality) bool {
        return self.reflexive and self.transitive and self.antisymmetric;
    }

    pub fn label(self: Causality) []const u8 {
        if (self.isPartialOrder()) return "porządek częściowy";
        if (self.transitive) return "praporządek (brak antysymetrii)";
        return "ani przechodnia, ani antysymetryczna";
    }
};

pub const Order = struct {
    f: DiagonalForm,
    /// Tolerancja rozpoznania wektora zerowego. W arytmetyce zmiennoprzecinkowej
    /// g(v,v) = 0 realizuje się jako |g| ~ eps, więc "zerowy" MUSI być
    /// pojęciem z tolerancją. W MRS tolerancja jest jawnym parametrem typu,
    /// a nie magiczną stałą rozsianą po kodzie.
    tol: f64 = 1e-9,
    /// Indeks wymiaru wybranego jako strzałka czasu. Przy p ≥ 2 to jest
    /// DODATKOWA struktura, której forma sama nie wyznacza.
    time_arrow: usize,

    pub fn init(s: Signature, time_arrow: usize) OrderError!Order {
        if (s.n() > MAX_DIM) return error.DimensionTooLarge;
        if (time_arrow >= s.n()) return error.TimeArrowOutOfRange;
        if (s.roles[time_arrow] != .temporal) return error.TimeArrowNotTemporal;
        return .{ .f = DiagonalForm.init(s), .time_arrow = time_arrow };
    }

    pub fn dim(self: Order) usize {
        return self.f.n();
    }

    pub fn signature(self: Order) Signature {
        return self.f.signature;
    }

    pub fn norm2(self: Order, v: []const f64) f64 {
        return self.f.eval(v);
    }

    /// Domknięty stożek przyszłości. Konwencja znaku wchodzi WYŁĄCZNIE
    /// przez `classify` — dzięki temu funkcja działa i dla (+,−,−,−),
    /// i dla (−,+,+,+).
    pub fn inFutureCone(self: Order, v: []const f64) bool {
        if (self.f.classifyTol(v, self.tol) == .spatial) return false;
        return v[self.time_arrow] >= 0.0;
    }

    /// Oriented relation: u ⪯ v.
    pub fn leq(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.inFutureCone(d[0..u.len]);
    }

    /// Separacja niezorientowana: g(v−u, v−u) ≤ 0.
    /// To NIE jest ta sama relacja co ⪯ — jest symetryczna na stożku
    /// zerowym, więc nie jest porządkiem.
    pub fn separation(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.f.classifyTol(d[0..u.len], self.tol) != .spatial;
    }

    /// Separacja zerowa: g(v−u, v−u) = 0 — relacja równoważności,
    /// której klasy to promienie świetlne.
    pub fn nullSeparated(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.f.classifyTol(d[0..u.len], self.tol) == .null_like;
    }

    // -- generatory wektorów w stożku ----------------------------------------

    /// Losowy wektor w stożku przyszłości. Konstrukcja jednolita dla obu
    /// konwencji znaku:
    ///   1. losujemy składowe poza strzałką czasu,
    ///   2. skalujemy grupę przestrzenną tak, by jej wkład do g wynosił
    ///      Σ_temporal_rest x² + 1 — wtedy wkład wszystkich wymiarów poza
    ///      strzałką jest równy dokładnie −s_arrow,
    ///   3. g = s_arrow(t² − 1), więc t = 1 daje wektor zerowy, a t = 2
    ///      wektor czasowy — w OBU konwencjach, bo "czasowy" znaczy
    ///      "składowa czasowa dominuje", a nie "g ma konkretny znak".
    pub fn randomConeVec(
        self: Order,
        out: *[MAX_DIM]f64,
        rnd: std.Random,
        mode: ConeMode,
    ) bool {
        const s = self.signature();
        const nn = s.n();
        const s_arrow = s.signAt(self.time_arrow);

        for (0..nn) |i| out[i] = (rnd.float(f64) - 0.5) * 2.0;
        out[self.time_arrow] = 0;

        var q_space: f64 = 0;
        var q_time_rest: f64 = 0;
        for (0..nn) |i| {
            if (i == self.time_arrow) continue;
            switch (s.roles[i]) {
                .spatial => q_space += out[i] * out[i],
                .temporal => q_time_rest += out[i] * out[i],
                .degenerate => {},
            }
        }

        if (q_space <= 1e-12) return false; // brak wymiaru przestrzennego = brak stożka
        const f = @sqrt((q_time_rest + 1.0) / q_space);
        for (0..nn) |i| {
            if (s.roles[i] == .spatial) out[i] *= f;
        }

        var other: f64 = 0;
        for (0..nn) |i| {
            if (i == self.time_arrow) continue;
            other += s.signAt(i) * out[i] * out[i];
        }
        if (@abs(other + s_arrow) > 1e-9) return false; // konstrukcja nie wyszła

        out[self.time_arrow] = switch (mode) {
            .null_future => 1.0,
            .timelike_future => 2.0,
        };

        const cls = self.f.classifyTol(out[0..nn], self.tol);
        return switch (mode) {
            .null_future => cls == .null_like,
            .timelike_future => cls == .temporal,
        };
    }

    /// Losowy wektor w stożku o losowym trybie (czasowy albo zerowy).
    /// Kontrprzykłady na przechodniość dla p ≥ 2 wymagają wektorów
    /// ZEROWYCH (dowód: suma dwóch czasowych o tym samym zwrocie
    /// strzałki zostaje w stożku), więc generator musi umieć oba.
    pub fn randomConeVecAny(self: Order, out: *[MAX_DIM]f64, rnd: std.Random) bool {
        const mode: ConeMode = if (rnd.float(f64) < 0.5) .timelike_future else .null_future;
        return self.randomConeVec(out, rnd, mode);
    }

    // -- wyszukiwanie świadków ----------------------------------------------

    /// Naruszenie przechodniości ⪯: u ⪯ w, w ⪯ v, ale u ⋠ v.
    pub fn searchTransitivityViolation(
        self: Order,
        trials: usize,
        rnd: std.Random,
    ) ?Witness {
        const nn = self.dim();
        var u: [MAX_DIM]f64 = undefined;
        var d1: [MAX_DIM]f64 = undefined;
        var d2: [MAX_DIM]f64 = undefined;
        var total: [MAX_DIM]f64 = undefined;

        for (0..trials) |_| {
            for (0..nn) |i| u[i] = (rnd.float(f64) - 0.5) * 4.0;
            if (!self.randomConeVecAny(&d1, rnd)) continue;
            if (!self.randomConeVecAny(&d2, rnd)) continue;
            for (0..nn) |i| total[i] = d1[i] + d2[i];
            if (!self.inFutureCone(total[0..nn])) {
                var wit = Witness{ .len = nn };
                @memcpy(wit.u[0..nn], u[0..nn]);
                for (0..nn) |i| {
                    wit.w[i] = u[i] + d1[i];
                    wit.v[i] = u[i] + d1[i] + d2[i];
                }
                return wit;
            }
        }
        return null;
    }

    /// Czy stożek przyszłości jest wypukły (próbkowanie).
    /// Wypukłość ⟺ przechodniość ⪯.
    pub fn coneIsConvex(self: Order, trials: usize, rnd: std.Random) bool {
        const nn = self.dim();
        var a: [MAX_DIM]f64 = undefined;
        var b: [MAX_DIM]f64 = undefined;
        var sum: [MAX_DIM]f64 = undefined;
        var checked: usize = 0;
        for (0..trials) |_| {
            if (!self.randomConeVecAny(&a, rnd)) continue;
            if (!self.randomConeVecAny(&b, rnd)) continue;
            for (0..nn) |i| sum[i] = a[i] + b[i];
            checked += 1;
            if (!self.inFutureCone(sum[0..nn])) return false;
        }
        return checked > 0;
    }

    /// Konstrukcyjny świadek nieprzechodniości dla p ≥ 2.
    /// Wymiary: t0 i t1 czasowe, s przestrzenny.
    ///   a   = e_t0 + e_s      → g = s_t − s_s = 0         (zerowy)
    ///   v−a = e_t1 + e_s      → g = 0                      (zerowy)
    ///   v   = e_t0 + e_t1 + 2·e_s → g = 2·s_t − 4·s_s < 0  (przestrzenny)
    /// Widać, że nierówność "2·s_t < 4·s_s" zachodzi w obu konwencjach,
    /// bo s_s = −s_t.
    pub fn canonicalWitness(s: Signature) OrderError!?Witness {
        if (s.p() < 2) return null;
        var idx_t: [2]usize = undefined;
        var nt: usize = 0;
        var idx_s: usize = std.math.maxInt(usize);
        for (0..s.n()) |i| {
            switch (s.roles[i]) {
                .temporal => {
                    if (nt < 2) {
                        idx_t[nt] = i;
                        nt += 1;
                    }
                },
                .spatial => if (idx_s == std.math.maxInt(usize)) {
                    idx_s = i;
                },
                .degenerate => {},
            }
        }
        if (nt < 2 or idx_s == std.math.maxInt(usize)) return null;

        var wit = Witness{ .len = s.n() };
        wit.w[idx_t[0]] = 1.0;
        wit.w[idx_s] = 1.0;
        wit.v[idx_t[0]] = 1.0;
        wit.v[idx_t[1]] = 1.0;
        wit.v[idx_s] = 2.0;
        return wit;
    }

    /// Świadek braku antysymetrii ⪯: u ≠ v z u ⪯ v i v ⪯ u.
    /// Istnieje dokładnie wtedy, gdy r > 0 albo p ≥ 2 (Twierdzenie 5.2).
    /// Zwraca `null` dla sygnatur lorentzowskich — i to jest zgodne
    /// z twierdzeniem, a nie przypadkiem.
    pub fn searchAntisymmetryViolation(s: Signature) OrderError!?Witness {
        const nn = s.n();
        var idx_t: [2]usize = undefined;
        var nt: usize = 0;
        var idx_s: usize = std.math.maxInt(usize);
        var idx_rad: usize = std.math.maxInt(usize);
        for (0..nn) |i| {
            switch (s.roles[i]) {
                .temporal => {
                    if (nt < 2) {
                        idx_t[nt] = i;
                        nt += 1;
                    }
                },
                .spatial => if (idx_s == std.math.maxInt(usize)) {
                    idx_s = i;
                },
                .degenerate => if (idx_rad == std.math.maxInt(usize)) {
                    idx_rad = i;
                },
            }
        }

        var wit = Witness{ .len = nn };
        if (idx_rad != std.math.maxInt(usize)) {
            // wektor jądra: niezerowy, zerowa norma, zerowa składowa strzałki
            wit.v[idx_rad] = 1.0;
            return wit;
        }
        if (nt >= 2 and idx_s != std.math.maxInt(usize)) {
            // d: strzałka 0, druga składowa czasowa 1, przestrzeń 0.5
            // → g = s_t·(1 − 0.25), klasa czasowa w obu konwencjach
            wit.v[idx_t[1]] = 1.0;
            wit.v[idx_s] = 0.5;
            return wit;
        }
        return null;
    }

    /// Pełna diagnoza: czy ⪯ jest porządkiem częściowym i co konkretnie
    /// zawodzi. Weryfikacja empiryczna Twierdzenia 5.2.
    pub fn diagnose(self: Order, trials: usize, rnd: std.Random) OrderError!Causality {
        const trans_viol = self.searchTransitivityViolation(trials, rnd) != null;
        const anti_viol = (try Order.searchAntisymmetryViolation(self.signature())) != null;
        return .{
            .reflexive = true,
            .transitive = !trans_viol,
            .antisymmetric = !anti_viol,
        };
    }
};

/// Eksperymentalna weryfikacja Twierdzenia 5.1.
pub const Transitivity = struct {
    violations: usize,
    trials: usize,
    checked: usize,
    convex: bool,
};

pub fn probeTransitivity(self: Order, trials: usize, rnd: std.Random) Transitivity {
    const nn = self.dim();
    var d1: [MAX_DIM]f64 = undefined;
    var d2: [MAX_DIM]f64 = undefined;
    var total: [MAX_DIM]f64 = undefined;
    var violations: usize = 0;
    var checked: usize = 0;
    for (0..trials) |_| {
        if (!self.randomConeVecAny(&d1, rnd)) continue;
        if (!self.randomConeVecAny(&d2, rnd)) continue;
        for (0..nn) |i| total[i] = d1[i] + d2[i];
        checked += 1;
        if (!self.inFutureCone(total[0..nn])) violations += 1;
    }
    return .{
        .violations = violations,
        .trials = trials,
        .checked = checked,
        .convex = violations == 0,
    };
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "Minkowski 3+1: ⪯ jest porządkiem częściowym" {
    var prng = std.Random.DefaultPrng.init(2024);
    const rnd = prng.random();
    const o = try Order.init(sig.minkowski_3_1, 0);

    const zero = [_]f64{ 0, 0, 0, 0 };
    try std.testing.expect(o.leq(&zero, &zero)); // zwrotność

    const res = probeTransitivity(o, 20_000, rnd);
    try std.testing.expectEqual(@as(usize, 0), res.violations);
    try std.testing.expect(res.checked > 15_000);
    try std.testing.expect(res.convex);

    // antysymetria: brak świadka, i to jest treść twierdzenia
    try std.testing.expectEqual(
        @as(?Witness, null),
        try Order.searchAntisymmetryViolation(sig.minkowski_3_1),
    );
    try std.testing.expectEqual(
        @as(?Witness, null),
        try Order.canonicalWitness(sig.minkowski_3_1),
    );

    const diag = try o.diagnose(20_000, rnd);
    try std.testing.expect(diag.isPartialOrder());
}

test "Twierdzenie 5.1: p = 2 łamie przechodniość" {
    const o = try Order.init(sig.two_times_2_1, 0);

    const wit = (try Order.canonicalWitness(sig.two_times_2_1)).?;
    const u = wit.uSlice();
    const wv = wit.wSlice();
    const v = wit.vSlice();

    try std.testing.expect(o.leq(u, wv)); // 0 ⪯ a
    try std.testing.expect(o.leq(wv, v)); // a ⪯ v
    try std.testing.expect(!o.leq(u, v)); // 0 ⋠ v  ← naruszenie

    // konkretne liczby: a zerowy, v−a zerowy, v przestrzenny
    try std.testing.expectEqual(form.Class.null_like, o.f.classify(wv));
    var diff: [3]f64 = undefined;
    for (0..3) |i| diff[i] = v[i] - wv[i];
    try std.testing.expectEqual(form.Class.null_like, o.f.classify(&diff));
    try std.testing.expectEqual(form.Class.spatial, o.f.classify(v));

    // niezależne potwierdzenie: losowe próby też znajdują naruszenie
    var prng = std.Random.DefaultPrng.init(31);
    try std.testing.expect(o.searchTransitivityViolation(5000, prng.random()) != null);
}

test "Twierdzenie 5.2: brak antysymetrii dla p ≥ 2 i dla r > 0" {
    // p = 2, r = 0
    {
        const o = try Order.init(sig.two_times_2_1, 0);
        const wit = (try Order.searchAntisymmetryViolation(sig.two_times_2_1)).?;
        const u = wit.uSlice();
        const v = wit.vSlice();
        try std.testing.expect(o.leq(u, v));
        try std.testing.expect(o.leq(v, u));
        try std.testing.expect(v[1] != 0.0); // u ≠ v
    }
    // r = 1
    {
        const o = try Order.init(sig.degenerate_2_1_1, 0);
        const wit = (try Order.searchAntisymmetryViolation(sig.degenerate_2_1_1)).?;
        const u = wit.uSlice();
        const v = wit.vSlice();
        try std.testing.expect(o.leq(u, v));
        try std.testing.expect(o.leq(v, u));
        try std.testing.expect(v[3] != 0.0); // wektor jądra, u ≠ v
    }
    // p = 1, r = 0: brak świadka — porządek częściowy
    {
        try std.testing.expectEqual(
            @as(?Witness, null),
            try Order.searchAntisymmetryViolation(sig.minkowski_3_1),
        );
    }
    // euklidesowa (p = 0): nie ma stożka, nie ma strzałki czasu
    try std.testing.expectError(
        error.TimeArrowNotTemporal,
        Order.init(sig.euclidean_4, 0),
    );
}

test "Twierdzenie 5.2 — pełna charakteryzacja: porządek ⟺ lorentzowska" {
    var prng = std.Random.DefaultPrng.init(4242);
    const cases = [_]sig.Signature{
        sig.minkowski_3_1,
        sig.minkowski_3_1_flipped,
        sig.minkowski_1_1,
        sig.two_times_2_1,
        sig.degenerate_2_1_1,
    };
    for (cases) |s| {
        const o = try Order.init(s, 0);
        const diag = try o.diagnose(4000, prng.random());
        try std.testing.expectEqual(s.isLorentzian(), diag.isPartialOrder());
    }
}

test "separacja zerowa jest symetryczna, więc nie jest porządkiem" {
    const o = try Order.init(sig.minkowski_1_1, 0);
    const u = [_]f64{ 0, 0 };
    const v = [_]f64{ 1, 1 }; // zerowy: g = 1 − 1 = 0

    // relacja niezorientowana: oba kierunki, u ≠ v → brak antysymetrii
    try std.testing.expect(o.separation(&u, &v));
    try std.testing.expect(o.separation(&v, &u));
    try std.testing.expect(o.nullSeparated(&u, &v));
    try std.testing.expect(u[0] != v[0]);

    // relacja zorientowana: v ⪯ u NIE zachodzi (strzałka czasu maleje)
    try std.testing.expect(o.leq(&u, &v));
    try std.testing.expect(!o.leq(&v, &u));

    // na ilorazie przez separację zerową porządek jest ostry
    const x = [_]f64{ 2, 0 };
    try std.testing.expect(o.leq(&u, &x));
    try std.testing.expect(!o.leq(&x, &u));
}

test "generator wektorów zerowych faktycznie daje g = 0" {
    var prng = std.Random.DefaultPrng.init(4242);
    const rnd = prng.random();
    const o = try Order.init(sig.minkowski_3_1, 0);
    var v: [MAX_DIM]f64 = undefined;
    var hits: usize = 0;
    for (0..1000) |_| {
        if (o.randomConeVec(&v, rnd, .null_future)) {
            try std.testing.expectApproxEqAbs(@as(f64, 0.0), o.norm2(v[0..4]), 1e-9);
            hits += 1;
        }
    }
    try std.testing.expect(hits > 900);
}

test "generator wektorów czasowych działa w obu konwencjach" {
    var prng = std.Random.DefaultPrng.init(555);
    const rnd = prng.random();
    const cases = [_]sig.Signature{
        sig.minkowski_3_1,
        sig.minkowski_3_1_flipped,
        sig.two_times_2_1,
    };
    for (cases) |s| {
        const o = try Order.init(s, 0);
        var v: [MAX_DIM]f64 = undefined;
        var hits: usize = 0;
        for (0..500) |_| {
            if (o.randomConeVec(&v, rnd, .timelike_future)) {
                try std.testing.expect(o.inFutureCone(v[0..s.n()]));
                try std.testing.expectEqual(form.Class.temporal, o.f.classify(v[0..s.n()]));
                hits += 1;
            }
        }
        try std.testing.expect(hits > 450);
    }
}

test "Twierdzenie 5.3: r > 0 degeneruje przyczynowość" {
    const o = try Order.init(sig.degenerate_2_1_1, 0);

    const rad = [_]f64{ 0, 0, 0, 1 }; // wektor z jądra
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), o.norm2(&rad), 1e-15);
    try std.testing.expect(o.nullSeparated(&rad, &[_]f64{ 0, 0, 0, 0 }));

    // przyczynowość przestaje odróżniać punkty: 0 ⪯ rad ⪯ 0, a rad ≠ 0
    const zero = [_]f64{ 0, 0, 0, 0 };
    try std.testing.expect(o.leq(&zero, &rad));
    try std.testing.expect(o.leq(&rad, &zero));

    var out: [4]f64 = undefined;
    try std.testing.expectError(error.DegenerateForm, o.f.raiseIndex(&rad, &out));
}

test "strzałka czasu musi być wymiarem czasowym" {
    try std.testing.expectError(
        error.TimeArrowNotTemporal,
        Order.init(sig.minkowski_3_1, 1),
    );
    try std.testing.expectError(
        error.TimeArrowOutOfRange,
        Order.init(sig.minkowski_3_1, 9),
    );
}
