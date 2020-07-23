# Load packages
using DifferentialEquations
using MAT
using DSP
using SparseArrays
using LinearAlgebra

# Set parameters
const taui = 1;
const taue = 1;
const aee = 1.0 / taue;
const aie = 1.4 / taui;
const aei = 0.7 / taue;
const aii = 0.7 / taui;
const aext = 0.9 / taue;
const kee = 1.0 * pi * taue;
const kie = 1.0 * pi * taue;
const kei = 2.0 * pi * taui;
const kii = 3.0 * pi * taui;
const vee = 10;
const vie = 8;
const vei = -8;
const vii = -12;
const de = 0.5;
const di = 0.5;
const kc = 25;
const ni = -50.0;
const ne = 25.0;

# Conduction velocity and offset delay
cv = 0.009;
od = 0;

# Load connectivity/perturbation/initial condition data
file = matopen("Cnorm.mat");
const C = kc * read(file, "C");
close(file);
file = matopen("pert.mat");
const pert = read(file, "pert");
close(file);
const pert=0.01*rand(1092);
file = matopen("h.mat");
const IC = read(file, "h");
close(file);

# Number of nodes
const Nn = size(C, 1);

# Node degrees
const Cn = append!([0], cumsum(sum(C .!= 0, dims = 2), dims = 1));

# Number of ODEs for each node
const Vn = cumsum(repeat([14], 1, Nn), dims = 2) .- 14;

h(p, t; idxs = nothing) = IC[idxs];

# Define functions
f(x, y, tau) =
    (1 / (pi * tau)) * (1 .- x .^ 2 .- y .^ 2) ./
    (1 .+ 2 * x .+ x .^ 2 .+ y .^ 2)
RF(x, y, n, d) =
    x .* y - y - d * (x .^ 2 - y .^ 2 + 2 * x .+ 1) / 2 + n * (-x .* y - y);
IF(x, y, n, d) =
    (-x .^ 2 + y .^ 2 + 2 * x .- 1 - d * (2 * x .* y + 2 * y) +
     n * (x .^ 2 - y .^ 2 + 2 * x .+ 1)) / 2;
RG(x, y, g, v) = sum((-x .* y - y) * v .* g - (x .^ 2 - y .^ 2 .- 1) .* g / 2);
IG(x, y, g, v) = sum((x .^ 2 - y .^ 2 + 2 * x .+ 1) * v .* g / 2 - x .* y .* g);

# Define (non-delay) ODE system
function CB_TMS(dydt, S, p, t)

    for n = 1:Nn
        Fn = zeros(14, 1)
        o = 14 * (n - 1)

        xe = S[o+1]
        ye = S[o+2]
        xi = S[o+3]
        yi = S[o+4]
        gii = S[o+5]
        gie = S[o+6]
        gei = S[o+7]
        gee = S[o+8]
        gext = S[o+9]
        hii = S[o+10]
        hie = S[o+11]
        hei = S[o+12]
        hee = S[o+13]
        hext = S[o+14]

        Fn[1] = (1 / taue) * (RF(xe, ye, ne, de) + RG(xe, ye, gei, vei) +
                 RG(xe, ye, gee, vee) + RG(xe, ye, gext, vee))
        Fn[2] = (1 / taue) * (IF(xe, ye, ne, de) + IG(xe, ye, gei, vei) +
                 IG(xe, ye, gee, vee) + IG(xe, ye, gext, vee))
        Fn[3] = (1 / taui) * (RF(xi, yi, ni, di) + RG(xi, yi, gii, vii) +
                 RG(xi, yi, gie, vie))
        Fn[4] = (1 / taui) * (IF(xi, yi, ni, di) + IG(xi, yi, gii, vii) +
                 IG(xi, yi, gie, vie))

        Fn[5:9] = [hii; hie; hei; hee; hext]
        Fn[10] = aii^2 * (kii * f(xi, yi, taui) - gii - 2 / aii * hii)
        Fn[11] = aie^2 * (kie * f(xe, ye, taue) - gie - 2 / aie * hie)
        Fn[12] = aei^2 * (kei * f(xi, yi, taui) - gei - 2 / aei * hei)
        Fn[13] = aee^2 * (kee * f(xe, ye, taue) - gee - 2 / aee * hee)
        Fn[14] = aext^2 *
                 (f(S[Vn.+1], S[Vn.+2], taue)*C[n:n, :]'.-gext.-2/aext*hext)[1]

        dydt[Vn[n]+1:Vn[n]+14] = Fn
    end
end

# Define DDE system
function CB_TMS_D(dydt, S, h, p, t)
    for n = 1:Nn
        Fn = zeros(14, 1)
        o = 14 * (n - 1)
        a = findall(!iszero, C[n, :])
        c = size(a)[1]

        xe = S[o+1]
        ye = S[o+2]
        xi = S[o+3]
        yi = S[o+4]
        gii = S[o+5]
        gie = S[o+6]
        gei = S[o+7]
        gee = S[o+8]
        gext = S[o+9]
        hii = S[o+10]
        hie = S[o+11]
        hei = S[o+12]
        hee = S[o+13]
        hext = S[o+14]

        Fn[1] = (1 / taue) * (RF(xe, ye, ne, de) + RG(xe, ye, gei, vei) +
                 RG(xe, ye, gee, vee) + RG(xe, ye, gext, vee))
        Fn[2] = (1 / taue) * (IF(xe, ye, ne, de) + IG(xe, ye, gei, vei) +
                 IG(xe, ye, gee, vee) + IG(xe, ye, gext, vee))
        Fn[3] = (1 / taui) * (RF(xi, yi, ni, di) + RG(xi, yi, gii, vii) +
                 RG(xi, yi, gie, vie))
        Fn[4] = (1 / taui) * (IF(xi, yi, ni, di) + IG(xi, yi, gii, vii) +
                 IG(xi, yi, gie, vie))

        Fn[5:9] = [hii; hie; hei; hee; hext]
        Fn[10] = aii^2 * (kii * f(xi, yi, taui) - gii - 2 / aii * hii)
        Fn[11] = aie^2 * (kie * f(xe, ye, taue) - gie - 2 / aie * hie)
        Fn[12] = aei^2 * (kei * f(xi, yi, taui) - gei - 2 / aei * hei)
        Fn[13] = aee^2 * (kee * f(xe, ye, taue) - gee - 2 / aee * hee)
        delay_f = zeros(c)
        for m = 1:c
            delay_f[m] = f(
                h(p, t - lags[Cn[n]+m]; idxs = Vn[a[m]] + 1),
                h(p, t - lags[Cn[n]+m]; idxs = Vn[a[m]] + 2),
                taue,
            )
        end

        Fn[14] = aext^2 * (sum(delay_f .* C[n:n, a]')-gext-2/aext*hext)[1]

        dydt[Vn[n]+1:Vn[n]+14] = Fn
    end
end

# Integrate ODE system (to steady state)
tspan = (0.0, 20.0);
u0 = IC;
prob = ODEProblem(CB_TMS, u0, tspan);
alg = MethodOfSteps(Tsit5());
sol = solve(
    prob,
    maxiters = 1e8,
    alg,
    progress = true,
    abstol = 1e-8,
    reltol = 1e-8,
);

file = matopen("steadystate.mat", "w");
write(file, "sol", sol.u);
close(file);

file = matopen("steadystate_time.mat", "w");
write(file, "t", sol.t);
close(file);

# Define delays
file = matopen("lags.mat");
lag = round.(read(file, "lags") / cv);
close(file);
const lags = lag + repeat([od], length(lag), 1);
lag = nothing;

# Integrate DDE system (with initial perturbation from steady state)
tspan = (0.0, 200.0);
u0 = sol.u[end]+pert;
prob = DDEProblem(CB_TMS_D, u0, h, tspan; constant_lags = lags);
sol = solve(
    prob,
    maxiters = 1e8,
    alg,
    progress = true,
    abstol = 1e-8,
    reltol = 1e-6,
    saveat=tspan[2]-15000:0.4:tspan[2],
    save_idxs = append!(
        round.(Int, Vn[1:end] .+ 1)[:, 1],
        round.(Int, Vn[1:end] .+ 2)[:, 1],
    ),
);

# Save time series of firing rate for excitatory populations
U = zeros(size(sol.u)[1], 2 * Nn);

for n = 1:size(sol.u)[1]
    U[n, :] = sol.u[n]
end

U = f(U[:, 1:Nn], U[:, Nn+1:end], taue);

file = matopen("timeseries.mat", "w");
write(file, "U", U);
close(file);

file = matopen("timeseries_t.mat", "w");
write(file, "t", sol.t);
close(file);
