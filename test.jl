using Revise
using FunctionTabulations
using FunctionTabulations: create_tabulation_1D_no_units
using Unitful
using SHA
using Test

#%%
c1 = code_lowered(f1)
c2 = code_lowered(f1)

h1 = bytes2hex(sha256(string(c1)))
h2 = bytes2hex(sha256(string(c2)))
h1==h2


code_lowered(f1)



@less create_tabulation_1D(f2, xmin=0.0, xmax=2pi, npoints=100)
hash("f1")
#%%
f1(x) = sin(x)^2 
f2(x) = sin(x)
x = collect(range(0.0,2pi, length=100))
itp = create_tabulation_1D(f1, x)
itp2 = create_tabulation_1D(f2, xmin=0.0, xmax=2pi, npoints=100)


isapprox(itp(1.0),  itp2(1.0))
#%%
f1b(x) = x^2
f2b(x) = x^2
x = collect(range(0.0,2.0, length=100))u"m"
itpb = create_tabulation_1D(f1b, x)
itp2b = create_tabulation_1D(f2b, xmin=0.0u"m", xmax=2.0u"m", npoints=100)


isapprox(itpb(1.0u"m"),  itp2b(1.0u"m"))
#%%

func_1d(x) = x
fw1 = FunctionTabulations.wrap_function_1D_remove_units(func_1d, u"m", u"m")
c1 = code_lowered(fw1) 
# create interpolations and test them

itp_1d_1 = create_tabulation_1D(
    func_1d,
    xmin = 0.0u"m",
    xmax = 3.0u"m",
    npoints = 100,
    x_scale = :linear,
    f_scale = :linear
)

itp_1d_1(1u"m")

func_1d(x) = x^2
fw2 = FunctionTabulations.wrap_function_1D_remove_units(func_1d, u"m", u"m^2")
c2 = code_lowered(fw2) 

c1 == c2

itp2 = create_tabulation_1D(
    func_1d,
    xmin = 0.0u"m",
    xmax = 3.0u"m",
    npoints = 100,
    x_scale = :linear,
    f_scale = :linear
)

itp2(2u"m")

@test_warn "The SHA for $func_1d did not match the one of the stored tabulated function. Please check if the function definition has changed" create_tabulation_1D(
    func_1d,
    xmin = 0.0u"m",
    xmax = 3.0u"m",
    npoints = 100,
    x_scale = :linear,
    f_scale = :linear
)
rm("func_1d_data.jld2")