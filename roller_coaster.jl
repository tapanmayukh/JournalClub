### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> title = "Roller Coaster Dynamics"
#> description = "From Newtonian Mechanics to Lagrangian Density: Analysing the roller coaster dynamics on different tracks."
#> date = "2026-09-30"
#> 
#>     [[frontmatter.author]]
#>     name = "Tapan Mayukh"
#>     [[frontmatter.author]]
#>     name = "Rohit Verma"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 25ee97bc-8efe-4471-9199-50291361d377
begin
    import Pkg;
    Pkg.activate(".");
end

# ╔═╡ 4466300d-4657-4b34-bd3c-272a19550f9b
using Plots, PlutoUI, PlotThemes, Revise, DifferentialEquations, Roots, LaTeXStrings;

# ╔═╡ ff3f2dee-a847-4df5-bf8c-6fc7c22ebe70
md"""
# Roller Coaster Dynamics
### From point particles to a continuum model using Lagrange density
"""

# ╔═╡ ab865418-5e3a-4062-ad56-e4df175f0128
html"""
<p align="right">&mdash; Presented By: <i>Rohit Verma</i> and <i>Tapan Mayukh</i></p>
"""

# ╔═╡ 8cbf05ba-b6a4-4e9c-ac38-b7a3963ca1b5


# ╔═╡ 1995dd8a-5ee0-4c26-8dfc-95135f36bf50


# ╔═╡ 09a2545f-7aa9-4a9c-a215-e1e710ccd99d


# ╔═╡ 55923503-d43f-48fb-b0e9-d44cea7043f2


# ╔═╡ 9bb46c1f-2b04-445b-a245-e2d12f14f94d
md"""
# Roller Coaster --- How does it work?
"""

# ╔═╡ 0a02ba54-5256-4b0e-9312-cabc6e8cce45
md"""

- Track is engineered to experience forces that you don't experience in daily life.

- Body induces stress response (adrenaline rush).

- At the end of ride, when your body feels safe, it releases Dopamine and Endorphins.
"""

# ╔═╡ d0151128-78e7-43a6-8b5c-e595e30dd25d


# ╔═╡ f738b9b4-ee60-4a2c-9cc2-a15beb7cbcc6
md"""
## Typical Roller Coaster trajectory elements
"""

# ╔═╡ 69bad9fd-c9c5-4699-93d1-103505692e63
track_list_md = md"""
There are many track types for a roller coaster:
1. Camelback
2. Loop-de-loop
3. In-line twist
4. Overbanked turn
5. Cobra roll
6. Top hat
7. More at [Coasterpedia - Roller Coaster Elements](https://coasterpedia.net/wiki/List_of_roller_coaster_elements)
""";

# ╔═╡ ba2126d7-43bc-4878-9f4a-fd6dfc2749b6
md"""
## Camelback Trajectory
"""

# ╔═╡ 517d9166-36e2-4fb6-b991-404b8023f811


# ╔═╡ 3a2f3ffc-189e-4450-966b-73173f8df444
md"""
## Camelback Trajectory --- Modelling
"""

# ╔═╡ a835de6d-e0ed-481a-a5d1-db0976e251c0


# ╔═╡ a1c7d6ac-bb2e-480c-8469-44bbb1c99e60
md"""
Show FWHM: $(@bind show_fwhm Switch())
"""

# ╔═╡ 1bf8a803-bb08-4e55-b014-991ec14c6824
md"""
## Single Wagon
"""

# ╔═╡ 1ce1a72c-85de-41b8-9de7-9c3548f94791


# ╔═╡ 6fa60588-963d-4b6b-95ff-37b5ba75b20c
md"""
## Single Wagon --- Minimum Climbing Velocity
"""

# ╔═╡ d5ce4638-3890-49ea-91f9-15b54c6e96e1


# ╔═╡ 31bb8c30-714c-4b72-a0fc-b9ae167a5fee
md"""
The gravitational force is conservative, and the Normal force is perpendicular to the motion, thus the total energy of the wagon is conserved. So, for the minimum velocity required to climb the hill:

```math
\frac{1}{2} m (\dot{x}_{\text{min}}^2 + \dot{z}_{\text{min}}^2) + m \, g \, z = m \, g \, H
```
"""

# ╔═╡ 12fe1903-e002-4749-be65-6554609320c7


# ╔═╡ cef6e64f-d9fd-4842-a7be-ff46f3787a8d
md"""
But ``\dot{z} = h' \, \dot{x}`` and ``z = h``, where ``h \equiv h(x)``, then:
```math
\dot{x}_{\text{min}} = \sqrt{\dfrac{2 \, g \, (H - h)}{1 + h'^2}}
```
"""

# ╔═╡ 63d25314-eee0-408f-b3a1-eefafa81a953


# ╔═╡ 51b15b08-75f2-4c22-b515-c0b1cb359a5f
md"""
As ``x \to \infty``, ``h \to 0`` and ``h' \to 0``, then:
```math
\dot{x}_{\text{min}} = \sqrt{2\, g \, H}
```
"""

# ╔═╡ ecf6d1bb-e076-4481-a8e1-59afeec5609f
md"""
## Single Wagon --- Equation of Motion
"""

# ╔═╡ af81083b-c29a-419b-90a7-3223a198fe87


# ╔═╡ 7a2df5f6-5ab9-4759-96ad-e60e44200fb1
md"""
**Force equations:**
```math
\begin{aligned}
m \, \ddot{x} &= - N \, \sin \theta \\[4pt]
m \, \ddot{z} &= N \, \cos \theta - m \, g
\end{aligned}
```
"""

# ╔═╡ a93a5477-2586-4b12-8841-77a1ba3a284c


# ╔═╡ bd3fa162-19b3-4ca8-84dc-4a7e662c0000
md"""
But,
```math
\begin{aligned}
\tan \theta &= \frac{dz}{dx} = h' \implies \dot{z} = h' \, \dot{x} \\[6pt]
\implies \ddot{z} &= h'' \, \dot{x}^2 + h' \, \ddot{x}
\end{aligned}
```
"""

# ╔═╡ c72f6a42-4215-4eee-a9c8-97cc62855384


# ╔═╡ 313af5c9-ee74-4655-ae4e-240ef9dc6a22
md"""
Thus,
```math
\begin{aligned}
\ddot{x} &= - h' \left(\dfrac{g + h'' \, \dot{x}^2}{1 + h'^2}\right) \\[6pt]
\ddot{z} &= \dfrac{h'' \, \dot{x}^2 - g \, h'^2}{1 + h'^2}
\end{aligned}
```
"""

# ╔═╡ 7302dd08-8601-451e-91f5-5411dd670750
md"""
## Single Wagon --- Plots
"""

# ╔═╡ ac303828-ea10-4266-a44a-826623740958
function plot_energy_single_wagon(sol, params; size_w=1000, size_h=700)
	kinetic_energy = 0.5 .* (sol[3, :].^2 + sol[4, :].^2)
	potential_energy = params.g .* sol[2, :]
	total_energy = kinetic_energy .+ potential_energy

	theme(:juno)
	plot([0, 0], [-10, 395]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol[1, :], kinetic_energy; lw=2, label="KE", color="#FE4365")
	plot!(sol[1, :], potential_energy; lw=2, label="PE", color="#eca25c")
	plot!(sol[1, :], total_energy; lw=2, label="Tot. E", color="#3f9778")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Energy Plot", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-10, 395)
end;

# ╔═╡ da9420e4-1ddc-400a-8993-b6fee5eda3ae
md"""
## Let's add another Wagon

(Both wagons are independent)
"""

# ╔═╡ 4163a1cf-8e8e-4337-9796-a9d7c6cf8608
function plot_separation_double_free_wagon(sol, params; size_w=1000, size_h=700)
	dx21 = sol[5, :] .- sol[1, :]
	dz21 = sol[6, :] .- sol[2, :]
	ds21 = sqrt.((dx21).^2 + (dz21).^2)

	theme(:juno)
	plot([0, 0], [-2, 2]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-0.5, 8], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol.t, ds21; lw=2, label=L"\mathbf{\Delta s_{2, 1}}", color="#FE4365")
	plot!(sol.t, dx21; linestyle=:dash, lw=2, label=L"\mathbf{\Delta x_{2, 1}}", color="#eca25c")
	plot!(sol.t, dz21; linestyle=:dash, lw=2, label=L"\mathbf{\Delta z_{2, 1}}", color="#3f9778")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Separation", titlelocation=:center)
	xlabel!("t")
	xlims!(-0.5, 8)
	ylims!(-2, 2)
end;

# ╔═╡ 79eb156c-7fd4-4588-acbc-cbdbeec95a18
md"""
## Two Wagons --- connected
"""

# ╔═╡ 0e78beae-b936-430c-8f08-11c6a48fb884
md"""
This gives us 3 constraint forces: ``N_{\text{front}}, N_{\text{rear}}, T``, which can be difficult to solve for using Newtonian Mechanics.
"""

# ╔═╡ afc900d4-557c-4d1e-91cd-ae89ff027291
md"""
## Two wagons --- fixed speed
(The Lagrangian way)
"""

# ╔═╡ 524258a3-440e-44f3-85b3-c28a58903836


# ╔═╡ cdccc6e1-c805-4d07-95fa-7b02bd51b7c9
md"""
Let's define the constraints first:
```math
g_i = z_i - h(x_i) = 0; \qquad \qquad g_3 = (\dot{x_1}^2 + \dot{z_1}^2) - (\dot{x_2}^2 + \dot{z_2}^2) = 0
```
"""

# ╔═╡ 39776d7a-1ba8-4821-9211-978e0ade0554
md"""
Lagrangian:
```math
ℒ = \dfrac{1}{2} m \, (\dot{x_1}^2 + \dot{z_1}^2 + \dot{x_2}^2 + \dot{z_2}^2) - m \, g \, (z_1 + z_2)
```
"""

# ╔═╡ 439eda1f-fdcb-4cb2-b06c-275ed6f1b16c
md"""
Lagrangian EOM:
```math
\dfrac{d}{dt}\dfrac{∂ℒ}{∂\dot{x_i}} - \dfrac{∂ ℒ}{∂x_i} = \sum_j \, λ_j \dfrac{∂ g_j}{∂ x_i}
```
"""

# ╔═╡ ad911bfc-3f35-4173-a849-765a0ee2c88d
L"""
\ddot{x_i} = - h_i' \left(\dfrac{g + 2 h_i'' \, \dot{x_i}^2}{2 (1 + h_i'^2)}\right) \qquad \qquad
\ddot{z_i} = \dfrac{2 h_i'' \, \dot{x_i}^2 - g \, h_i'^2}{2(1 + h_i'^2)}
"""

# ╔═╡ 85f5bb91-2b5c-4a36-a711-738264410ae0
md"""
## Two wagons --- fixed speed
"""

# ╔═╡ 1d4d6278-8e6b-4d93-a0eb-586202d4e4c1
function plot_separation_double_taut_wagon(sol, params; size_w=1000, size_h=700)
	dx21 = sol[5, :] .- sol[1, :]
	dz21 = sol[6, :] .- sol[2, :]
	ds21 = sqrt.((dx21).^2 + (dz21).^2)

	sep_max = ceil(maximum(ds21))

	theme(:juno)
	plot([0, 0], [-sep_max, sep_max]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-0.5, 8], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol.t, ds21; lw=2, label=L"\mathbf{\Delta s_{21}}", color="#FE4365")
	plot!(sol.t, dx21; linestyle=:dash, lw=2, label=L"\mathbf{\Delta x_{21}}", color="#eca25c")
	plot!(sol.t, dz21; linestyle=:dash, lw=2, label=L"\mathbf{\Delta z_{21}}", color="#3f9778")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	title!("Separation", titlelocation=:center)
	xlabel!("t")
	xlims!(-0.5, 8)
	ylims!(-sep_max, sep_max)
end;

# ╔═╡ 149b825b-3ea7-4aeb-9d9e-94291ddd0bf7
md"""
## Two Wagons --- connected by a spring
"""

# ╔═╡ 9d6624a4-3142-41bb-aebf-1ebbf66a20a7


# ╔═╡ 1ce2148a-49a8-4fb5-822d-e4bf4b53fdb7
md"""
Let's define the constraints first:
```math
g_i = z_i - h(x_i) = 0
```
"""

# ╔═╡ 36fd7fe2-e34b-4118-9bb6-a797ec963bfc
md"""
Lagrangian:
```math
\begin{aligned}
ℒ &= \dfrac{1}{2} m \, \sum_{i=1}^N \, (\dot{x_i}^2 + \dot{z_i}^2) - m \,g \, \sum_{i=1}^N \, z_i - \dfrac{1}{2} \sum_{i=1}^{N-1} \, k \left(\Delta_{i+1, i} - s_0\right) \\[4pt]
\Delta_{i+1, i} &= \sqrt{(x_{i+1} - x_{i})^2 + (z_{i+1}- z_{i})^2}
\end{aligned}
```
"""

# ╔═╡ 72dfbb63-87c6-4218-9237-0a1b8e56493a
md"""
Then:
```math
\begin{aligned}
\ddot{x}_2 &= \dfrac{k S_{21} (x_2 - x_1) - h_1' (h_1'' \dot{x}_1^2 + g - k S_{21} (z_2 - z_1)}{1 + h_1'^2}; \qquad S_{21} = \dfrac{\Delta_{21} - s_0}{\Delta_{21}} \\[4pt]
\ddot{x}_2 &= \dfrac{- k S_{21} (x_2 - x_1) - h_2' (h_2'' \dot{x}_2^2 + g + k S_{21} (z_2 - z_1))}{1 + h_2'^2}
\end{aligned}
```
"""

# ╔═╡ 730c2c53-bc09-4783-a40a-2bdc5c21c25c
md"""
## Two Wagons --- connected by a spring
(a hard spring such that ``\frac{k}{m} \geq 10 \frac{g}{s_0}``)
"""

# ╔═╡ 158ad312-bbf0-492b-aede-23ae67ee0198
md"""
## Two Wagons --- connected by a spring
(a hard spring such that ``\frac{k}{m} \geq 10 \frac{g}{s_0}``)
"""

# ╔═╡ 2b82e765-6017-4d92-8c0d-b4aab70a9110


# ╔═╡ 052494d1-526f-45ef-a27c-f2e7029d5b4d
function plot_velocity_double_spring_wagon(sol, params; size_w=1000, size_h=700, time=true)
	v1 = sqrt.(sol[3, :].^2 .+ sol[4, :].^2)
	v2 = sqrt.(sol[7, :].^2 .+ sol[8, :].^2)
	
	v1_max = ceil(maximum(v1))
	v2_max = ceil(maximum(v2))

	if v1_max > v2_max
		v_max = v1_max
	else
		v_max = v2_max
	end
	
	theme(:juno)
	plot([0, 0], [-5, v_max + 5]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	
	if time
		plot!([-0.5, 8], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
		plot!(sol.t, v1; lw=2, label=L"\mathbf{v_1}", color="#FE4365")
		plot!(sol.t, v2; lw=2, label=L"\mathbf{v_2}", color="#eca25c")
		xlabel!("t")
		xlims!(-0.5, 8)
		ylims!(-5, v_max + 5)
	else
		plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
		plot!(sol[1, :], v1; lw=2, label=L"\mathbf{v_1}", color="#FE4365")
		plot!(sol[5, :], v2; lw=2, label=L"\mathbf{v_2}", color="#eca25c")
		xlabel!("x")
		xlims!(-50, 50)
		ylims!(-5, v_max + 5)
	end

	title!("Velocity", titlelocation=:center)
end;

# ╔═╡ 6ca184b5-133b-4c9a-a691-3a59da8663f5
function plot_separation_double_spring_wagon(sol, params; size_w=1000, size_h=700)
	dx21 = sol[5, :] .- sol[1, :]
	dz21 = sol[6, :] .- sol[2, :]
	ds21 = sqrt.((dx21).^2 + (dz21).^2)

	sep_max = ceil(maximum(ds21))

	theme(:juno)
	plot([0, 0], [0.95, 1.05]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-0.5, 8], [1, 1]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol.t, ds21; lw=2, label=L"\mathbf{\Delta s_{21}}", color="#FE4365")
	# plot!(sol.t, dx21; linestyle=:dash, lw=2, label="\\Delta x_{21}", color="#eca25c")
	# plot!(sol.t, dz21; linestyle=:dash, lw=2, label="\\Delta z_{21}", color="#3f9778")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	title!("Separation", titlelocation=:center)
	xlabel!("t")
	xlims!(-0.5, 8)
	ylims!(0.95, 1.05)
end;

# ╔═╡ f195aae0-6e8b-41b8-a32b-4325bbf2f41b
md"""
## N wagons --- connected by springs
"""

# ╔═╡ 257585f4-9231-49ab-a787-ae8fa01ceb6a
function plot_velocity_n_spring_wagon(sol, params; size_w=1000, size_h=700, time=true)
	n = params.n
	middle = 4 * floor(Int, n/2)
	
	v1 = sqrt.(sol[3, :].^2 .+ sol[4, :].^2)
	v2 = sqrt.(sol[middle+3, :].^2 .+ sol[middle+4, :].^2)
	v3 = sqrt.(sol[end-1, :].^2 .+ sol[end, :].^2)
	
	v1_max = ceil(maximum(v1))
	v2_max = ceil(maximum(v2))
	v3_max = ceil(maximum(v3))

	v_max = maximum([v1_max, v2_max, v3_max])
	
	theme(:juno)
	plot([0, 0], [-5, v_max + 5]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	
	if time
		plot!([-0.5, 8], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
		plot!(sol.t, v1; lw=2, label=L"v_{front}", color="#FE4365")
		plot!(sol.t, v2; lw=2, label=L"v_{middle}", color="#eca25c")
		plot!(sol.t, v3; lw=2, label=L"v_{rear}", color="#3f9778")
		xlabel!("t")
		xlims!(-0.5, 8)
		ylims!(-5, v_max + 5)
	else
		plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
		plot!(sol[1, :], v1; lw=2, label=L"v_{front}", color="#FE4365")
		plot!(sol[middle+1, :], v2; lw=2, label=L"v_{middle}", color="#eca25c")
		plot!(sol[end-3, :], v2; lw=2, label=L"v_{rear}", color="#3f9778")
		xlabel!("x")
		xlims!(-50, 50)
		ylims!(-5, v_max + 5)
	end

	title!("Velocity", titlelocation=:center)
end;

# ╔═╡ d93eeb22-e1bb-478c-a97c-7327f558800c
function plot_separation_n_spring_wagon(sol, params; size_w=1000, size_h=700, fixed=false)
	n = params.n
	s0 = params.sep

	if fixed
		s0 /= (n - 1)
	end
	
	dx21 = sol[5, :] .- sol[1, :]
	dz21 = sol[6, :] .- sol[2, :]
	ds21 = sqrt.((dx21).^2 + (dz21).^2)

	dx32 = sol[end-3, :] .- sol[end-7, :]
	dz32 = sol[end-2, :] .- sol[end-6, :]
	ds32 = sqrt.((dx32).^2 + (dz32).^2)

	sep_max = maximum([ceil(maximum(ds21)), ceil(maximum(ds32))])
	
	theme(:juno)
	plot([0, 0], [s0-0.1*s0, s0+0.1*s0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-0.5, 8], [s0, s0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol.t, ds21; lw=2, label="wagon2 - wagon1", color="#FE4365")
	plot!(sol.t, ds32; lw=2, label="wagon$(n) - wagon$(n-1)", color="#eca25c")
	# plot!(sol.t, dx21; linestyle=:dash, lw=2, label="\\Delta x_{21}", color="#eca25c")
	# plot!(sol.t, dz21; linestyle=:dash, lw=2, label="\\Delta z_{21}", color="#3f9778")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	title!("Separation", titlelocation=:center)
	xlabel!("t")
	xlims!(-0.5, 8)
	ylims!(s0-0.1*s0, s0+0.1*s0)
end;

# ╔═╡ e43f9ff7-7074-4129-a579-121fccd2f570
md"""
## N wagons --- connected by springs (Fixed Initial Length)
"""

# ╔═╡ ef2b1405-63c3-4f52-aaef-4adf9ae37ddb
md"""
## Wagon Continuum --- Lagrangian
"""

# ╔═╡ 54d48c4f-5624-47e4-b752-fd161857fac0


# ╔═╡ 95f05467-2968-4fb5-b79d-353e11b5eaa9
md"""
In continuum limit, ``s_0 = dl \to 0`` and ``N \to \infty``. Then, with ``k \, s_0 = \Phi`` and ``m = \rho s_0``:

```math
T = \dfrac{1}{2} \rho \sum_i (\dot{x_i}^2 + \dot{z_i}^2) \, s_0 = \dfrac{1}{2} \int_0^L  \rho \,\left(1 + h'^2\right) \left(\dfrac{\partial x}{\partial t}\right)^2 dl
```
"""

# ╔═╡ 2ccbfaf8-6283-49fd-9abc-f1127eddb11f


# ╔═╡ 1f7a9428-8c7b-48c1-928c-b4f33909a79a
md"""
Also:

```math
\begin{aligned}
U &= \rho g \sum_i z_i s_0 + \dfrac{1}{2} \dfrac{\Phi}{s_0} \sum_i \left(\sqrt{(x_{i+1} - x_i)^2 + (z_{i+1} - z_i)^2} - s_0\right)^2 \\[4pt]
&= \int_0^L \rho \, g \, h \, dl + \dfrac{1}{2} \Phi \, \int_0^L \left(\sqrt{1 + h'^2} \, \dfrac{∂x}{∂l} - 1\right)^2 dl
\end{aligned}
```
"""

# ╔═╡ aa9be589-48c6-41ee-bb29-39114d6582c5


# ╔═╡ 88426bfc-4043-43c2-92a0-44f7a22fee4e
md"""
Thus:

```math
\int_0^L dl \, ℒ = \int_0^L dl \, ℒ\left(x, \dfrac{∂x}{∂l}, \dfrac{∂x}{∂t}\right) = T - U
```
"""

# ╔═╡ ff3c39a4-a447-4b34-97f1-1ff2a7736e69
md"""
## Wagon Continuum --- Action
"""

# ╔═╡ cb9f10ce-f0d2-49da-b835-8e663545c116


# ╔═╡ e750ee9e-cce0-42fd-9e13-a01750d9fbde
md"""
```math
A = \int_{t_1}^{t_2} dt \int_0^L dl \, ℒ\left(x, \dfrac{∂x}{∂l}, \dfrac{∂x}{∂t}\right)
```
"""

# ╔═╡ 142ba9b5-9239-40ac-8d35-bc1c558d3016
md"""
Then:
```math
\delta A = \int_{t_1}^{t_2} \int_0^L dl \left[ℒ\left(x + \delta x, \dfrac{∂x}{∂l} + \dfrac{∂\, \delta x}{∂l}, \dfrac{∂x}{∂t} + \dfrac{∂\, \delta x}{∂t}\right) - ℒ\left(x, \dfrac{∂x}{∂l}, \dfrac{∂x}{∂t}\right) \right]
```
"""

# ╔═╡ 6bf2a94d-ff43-41b4-a867-c148cf8abbd6
md"""
Equation of motion:
```math
\dfrac{\partial}{\partial t} \dfrac{\partial ℒ}{\partial \, (\partial_t x)} + \dfrac{\partial}{\partial l} \dfrac{\partial ℒ}{\partial \, (\partial_l x)} - \dfrac{\partial ℒ}{\partial x} = 0
```
"""

# ╔═╡ 7acb55f0-1632-446f-b173-f15e6f56f67a
md"""
Thus:
```math
\dfrac{∂^2 x}{\partial t^2} = - \dfrac{h' h''}{1 + h'^2} \left(\dfrac{∂ x}{\partial t}\right)^2 - \dfrac{g h'}{1 + h'^2} + \dfrac{\Phi}{\rho} \dfrac{h' h''}{1 + h'^2} \left(\dfrac{∂ x}{\partial l}\right)^2 + \dfrac{\Phi}{\rho} \left(\dfrac{∂^2 x}{∂l^2}\right)
```
"""

# ╔═╡ c82157cd-b834-4871-8cfa-d5f32af7182c
md"""
## Wagon Continuum --- Solved
"""

# ╔═╡ efd5daf8-5040-4c3e-96a9-ea71268ec43a


# ╔═╡ f41a230e-8425-441c-a701-647fd3ede72f
md"""
%trajectory: $(@bind percent_trajectory Slider(0.001:0.001:0.55; default=0.05, show_value=true))
"""

# ╔═╡ 17b42861-dd59-4a3c-a6a3-6ebe6806a596
md"""
## Loop-de-Loop
Also known as a Maclaurin Trisectrix
"""

# ╔═╡ f5732129-f37e-4300-b31d-8d1ec72ea697


# ╔═╡ 67bbab74-fb4e-46a8-8ed4-12e8ef298059
md"""
## Camelback Trajectory --- Modelling
"""

# ╔═╡ 7764737a-d571-46b2-9617-d89d0d10cb04


# ╔═╡ 4280d95b-88a5-4a4e-822a-fe8fc4e2b0b0
md"""
## Single Wagon --- Equation of Motion
"""

# ╔═╡ 0ed9822d-3945-4525-adea-6042c2c4acfd
md"""
**Force equations:**
```math
\begin{aligned}
m \, \ddot{x} &= - N \, \sin \theta \\[4pt]
m \, \ddot{z} &= N \, \cos \theta - m \, g
\end{aligned}
```
"""

# ╔═╡ 8625a145-6ceb-43fb-bdb8-1368e0ad62e7
md"""
But,
```math
\begin{aligned}
\tan \theta &= \frac{dz}{dx} = \frac{dz/du}{dx/du} = \frac{z'}{x'} = h' \\[6pt]
\dot{z} &= z' \, \dot{u} \\[6pt]
\dot{x} &= x' \, \dot{u} \\[6pt]
\ddot{z} &= z'' \, \dot{u}^2 + z' \, \ddot{u} \\[6pt]
\ddot{x} &= x'' \, \dot{u}^2 + x' \, \ddot{u}
\end{aligned}
```
"""

# ╔═╡ 67ea193d-ace4-4264-b20f-35724ff5db46
md"""
Here,
```math
\ddot{u} = - \dfrac{\dot{u}^2 (x' \, x'' + z' \, z'') + g \, z'}{x'^2 + z'^2}
```
"""

# ╔═╡ 30db840c-6393-4922-82c2-a147c09b24e4
md"""
## Single Wagon --- Plots
"""

# ╔═╡ 572cf4cd-e2ef-4469-a86d-bb757048f527
begin
	trisec_x(u, b) = b * u * (u^2 - 3) / (u^2 + 1)
	trisec_z(u, b) = - b * (u^2 - 3) / (u^2 + 1)

	trisec_dx(u, b) = b * (u^4 + 6u^2 - 3) / (u^2 + 1)^2
	trisec_d2x(u, b) = 8b * u * (3 - u^2) / (u^2 + 1)^3

	trisec_dz(u, b) = - 8b * u / (u^2 + 1)^2
	trisec_d2z(u, b) = 8b * (3u^2 - 1) / (u^2 + 1)^3
end;

# ╔═╡ 9ad4900b-c1c8-4fab-8de8-005caa983b53
function single_wagon_trisectrix_update!(du, u, p, t)
	u1, du1 = u

	dx = trisec_dx(u1, p.b)
	d2x = trisec_d2x(u1, p.b)
	dz = trisec_dz(u1, p.b)
	d2z = trisec_d2z(u1, p.b)

	du[1] = du1
	du[2] = -(du1^2 * (dx*d2x + dz*d2z) + p.g * dz) / (dx^2 + dz^2)
end;

# ╔═╡ cd4addf0-2a3a-4aa2-82ad-f17f726e143d
function plot_normal_single_wagon_trisec(sol, params; size_w=1000, size_h=700)
	dx = trisec_dx.(sol[1, :], params.b)
	dz = trisec_dz.(sol[1, :], params.b)

	d2x = trisec_d2x.(sol[1, :], params.b)
	d2z = trisec_d2z.(sol[1, :], params.b)

	x = trisec_x.(sol[1, :], params.b)
	z = trisec_z.(sol[1, :], params.b)
	
	normal = sol[2, :].^2 .* (dx .* d2z .- dz .* d2x)
	normal ./= sqrt.(dx.^2 .+ dz.^2)
	
	theme(:juno)
	plot([0, 0], [-10, 40.0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(x, normal; lw=2, label="Normal", color="#FE4365")
	plot!(x, z; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Normal Force", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-10, 40)
end;

# ╔═╡ 24685f7d-b140-44a9-b9cf-d6d9a9471142
function plot_energy_single_wagon_trisec(sol, params; size_w=1000, size_h=700)
	dx = trisec_dx.(sol[1, :], params.b)
	dz = trisec_dz.(sol[1, :], params.b)

	x = trisec_x.(sol[1, :], params.b)
	z = trisec_z.(sol[1, :], params.b)
	
	kinetic_energy = 0.5 .* (dx.^2 + dz.^2) .* (sol[2, :]).^2
	potential_energy = params.g .* (z .+ params.b)
	total_energy = kinetic_energy .+ potential_energy

	theme(:juno)
	plot([0, 0], [-10, 600]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(x, kinetic_energy; lw=2, label="KE", color="#FE4365")
	plot!(x, potential_energy; lw=2, label="PE", color="#eca25c")
	plot!(x, total_energy; lw=2, label="Tot. E", color="#3f9778")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Energy Plot", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-10, 600)
end;

# ╔═╡ faf9e517-1f09-4ea9-8e6e-e5f13863e1a9
md"""
# Conclusion

So, what seats are the best in a roller coaster?

Middle one, of course. It's the safest, and better safe than sorry.
"""

# ╔═╡ 58e1580d-b594-4c90-aa44-897399e323f0


# ╔═╡ e753b964-95a0-4614-830d-0b1b4ebef5c6
md"""
# Thank You
"""

# ╔═╡ 22a4e90d-17c2-4187-8d8a-aac8f3601975
md"""
# Appendix
"""

# ╔═╡ 61e4180b-f082-4174-b06b-6ed2d1759f4d
md"""
## Presentation Helpers
"""

# ╔═╡ 05f46635-15bd-4132-b6cb-cee94e6f2ee3
import PlutoCSSMacros;

# ╔═╡ 81b8bfb2-5441-4e16-9e04-739c44c72a71
@PlutoCSSMacros.css_font_small md"""[^Ref]: Michael Kaschke, Holger Cartarius; Roller coaster dynamics—From point particles to a continuum model using Lagrange density. Am. J. Phys. 1 June 2026; 94 (6): 453–459. [https://doi.org/10.1119/5.0237959](https://doi.org/10.1119/5.0237959)
"""

# ╔═╡ 71b1e8fc-e5dc-41a1-8220-a36036bfd3ba
coaster_img_html = PlutoCSSMacros.css_vertical(
	Resource("https://upload.wikimedia.org/wikipedia/commons/2/2e/BlueHawk_SFOG_1.jpg", :height => 400, :alt => "Blue Hawk, GA"),
	md"*Blue Hawk, GA* - Wikimedia"
);

# ╔═╡ 851c6405-c9de-49a1-9caa-1f2c66423012
PlutoCSSMacros.css_horizontal(track_list_md, coaster_img_html; width=[40, 60])

# ╔═╡ fd72801e-3131-4950-8185-bc382d1e4618
PlutoCSSMacros.css_vertical(
	Resource("https://upload.wikimedia.org/wikipedia/commons/a/a3/Europa-Park_-_Silver_Star_%2836%29.JPG", :height => 600, :alt => "Silver Star, BaWü"),
	md"*Silver Star, BaWü* - Wikimedia"
)

# ╔═╡ 735e1043-20c8-4e46-a576-80eab7c8707e
track_gaussian_html = PlutoCSSMacros.css_vertical(
	Resource("https://upload.wikimedia.org/wikipedia/commons/a/a3/Europa-Park_-_Silver_Star_%2836%29.JPG", :width => "100%", :alt => "Silver Star, BaWü"),
	md"*Silver Star, BaWü* - Wikimedia"
);

# ╔═╡ 5441a648-a176-4dc1-8f58-dd24e6b459c6
PlutoCSSMacros.css_vertical(
	Resource("https://upload.wikimedia.org/wikipedia/commons/f/f4/Shockwave_coaster_sfot.jpg?utm_source=en.wikipedia.org&utm_campaign=imageinfo&utm_content=original", :height => 600, :alt => "Shockwave, TX"),
	md"*Shockwave, TX* - Wikimedia"
)

# ╔═╡ cf860bf9-fca1-4fd3-a1f8-465cf19f0273
track_trisectrix_html = PlutoCSSMacros.css_vertical(
	Resource("https://upload.wikimedia.org/wikipedia/commons/f/f4/Shockwave_coaster_sfot.jpg?utm_source=en.wikipedia.org&utm_campaign=imageinfo&utm_content=original", :width => "100%", :alt => "Silver Star, BaWü"),
	md"*Shockwave, TX* - Wikimedia"
);

# ╔═╡ 662be44c-48d9-4cd6-9993-984e67e74869
PlutoCSSMacros.add_global_css()

# ╔═╡ 2e65a895-6053-4768-9ba3-3115fef2ddac
PlutoCSSMacros.set_colors()

# ╔═╡ 021eb2b3-6243-492c-8140-702a05d05350
md"""
## Track Functions
"""

# ╔═╡ 1dd9a5c4-70c4-4435-a679-9f78716a8c7f
function α(FWHM)
	2.0 * sqrt(log(2.0)) / FWHM
end;

# ╔═╡ ae6d3d87-97ba-4d32-a560-004a94c6d506
begin
	FWHM = 20.0
	H = 30.0
	g = 9.81

	params = (
		FWHM = FWHM,
		H = H,
		g = g
	)

	alpha = α(params.FWHM)
	x0 = -4.0 / alpha
end;

# ╔═╡ d45ebb00-2292-4c29-a13f-29f080774585
function gaussian_h(x; FWHM=20.0, H=30.0)
	alpha = α(FWHM)
	H .* exp.(-alpha^2 .* x.^2)
end;

# ╔═╡ f1d5822e-d077-40d1-b910-7098fd133a7d
function plot_trajectory_continuum(sol_t, sol, params; size_w=1000, size_h=700, percent_t=0.1)
	n = params.n
	s0 = params.sep / (n - 1)

	x1 = sol[1, :]
	z1 = gaussian_h(sol[1, :]; FWHM=params.FWHM, H=params.H)

	nt = length(sol_t)
	t_curr = round(Int, percent_t * nt)
	x_curr = sol[:, t_curr]
	z_curr = gaussian_h(sol[:, t_curr]; FWHM=params.FWHM, H=params.H)
	
	theme(:juno)
	plot([0, 0], [-10, 50]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	# plot!(sol_t, ds21; lw=2, label="wagon2 - wagon1", color="#FE4365")
	# plot!(sol_t, ds32; lw=2, label="wagon$(n) - wagon$(n-1)", color="#eca25c")
	plot!(x1, z1; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(sol[:, t_curr], z_curr; lw=5, label=nothing, color="#FE4365")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	title!("Trajectory", titlelocation=:center)
	xlabel!("x")
	xlims!(-50, 50)
	ylims!(-10, 50)
end;

# ╔═╡ 10884f5a-1af7-4448-9d25-3efc22aff04e
function gaussian_h_prime(x; FWHM=20.0, H=30.0)
	alpha = α(FWHM)
	-2.0 * alpha^2 .* x .* gaussian_h.(x)
end;

# ╔═╡ 80e4ba7d-9493-41c0-a177-b06c988a86bc
function gaussian_h_double_prime(x; FWHM=20.0, H=30.0)
	alpha = α(FWHM)
	gaussian_h.(x) .* (4.0 * alpha^4 .* x.^2 .- 2.0 * alpha^2)
end;

# ╔═╡ 69f15733-5946-40ef-9c37-15b5f4c6b73b
function single_wagon_update!(du, u, p, t)
	x, z, vx, vz = u

	f₁ = gaussian_h_prime(x; FWHM=p.FWHM, H=p.H)
	f₂ = gaussian_h_double_prime(x; FWHM=p.FWHM, H=p.H)

	du[1] = vx
	du[2] = vz
	du[3] = - f₁ * (p.g + f₂ * vx^2) / (1 + f₁^2)
	du[4] = (f₂ * vx^2 - p.g * f₁^2) / (1 + f₁^2)
end;

# ╔═╡ b7b7ff06-32e9-45f4-b043-b32a1c6f3f01
function plot_normal_single_wagon(sol, params; size_w=1000, size_h=700)
	normal = params.g .+ gaussian_h_double_prime(sol[1, :]; FWHM=params.FWHM, H=params.H) .* sol[3, :].^2
	normal ./= sqrt.(1 .+ gaussian_h_prime(sol[1, :]; FWHM=params.FWHM, H=params.H).^2)
	
	theme(:juno)
	plot([0, 0], [-10, 90.0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol[1, :], normal; lw=2, label="Normal", color="#FE4365")
	plot!(sol[1, :], sol[2, :]; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Normal Force", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-10, 90)
end;

# ╔═╡ ddbf78aa-29ce-41e2-ba61-45aab0d9c942
function plot_normal_double_free_wagon(sol, params; size_w=1000, size_h=700)
	normal1 = params.g .+ gaussian_h_double_prime(sol[1, :]; FWHM=params.FWHM, H=params.H) .* sol[3, :].^2
	normal1 ./= sqrt.(1 .+ gaussian_h_prime(sol[1, :]; FWHM=params.FWHM, H=params.H).^2)
	normal2 = params.g .+ gaussian_h_double_prime(sol[5, :]; FWHM=params.FWHM, H=params.H) .* sol[7, :].^2
	normal2 ./= sqrt.(1 .+ gaussian_h_prime(sol[5, :]; FWHM=params.FWHM, H=params.H).^2)
	
	theme(:juno)
	plot([0, 0], [-10, 90.0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol[1, :], normal1; lw=2, label="Normal - Front Wagon", color="#FE4365")
	plot!(sol[5, :], normal2; lw=2, label="Normal - Rear Wagon", color="#eca25c")
	plot!(sol[1, :], sol[2, :]; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Normal Force", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-10, 90)
end;

# ╔═╡ 8b2b545c-5ed0-44be-a881-63096994ba62
function double_free_wagon_update!(du, u, p, t)
	x1, z1, vx1, vz1, x2, z2, vx2, vz2 = u

	x1f₁ = gaussian_h_prime(x1; FWHM=p.FWHM, H=p.H)
	x1f₂ = gaussian_h_double_prime(x1; FWHM=p.FWHM, H=p.H)

	x2f₁ = gaussian_h_prime(x2; FWHM=p.FWHM, H=p.H)
	x2f₂ = gaussian_h_double_prime(x2; FWHM=p.FWHM, H=p.H)

	du[1] = vx1
	du[2] = vz1
	du[3] = - x1f₁ * (p.g + x1f₂ * vx1^2) / (1 + x1f₁^2)
	du[4] = (x1f₂ * vx1^2 - p.g * x1f₁^2) / (1 + x1f₁^2)

	du[5] = vx2
	du[6] = vz2
	du[7] = - x2f₁ * (p.g + x2f₂ * vx2^2) / (1 + x2f₁^2)
	du[8] = (x2f₂ * vx2^2 - p.g * x2f₁^2) / (1 + x2f₁^2)
end;

# ╔═╡ 078dbc7f-e3c4-44fa-aa86-4de02574a0a5
function double_taut_wagon_update!(du, u, p, t)
	x1, z1, vx1, vz1, x2, z2, vx2, vz2 = u

	x1f₁ = gaussian_h_prime(x1; FWHM=p.FWHM, H=p.H)
	x1f₂ = gaussian_h_double_prime(x1; FWHM=p.FWHM, H=p.H)

	x2f₁ = gaussian_h_prime(x2; FWHM=p.FWHM, H=p.H)
	x2f₂ = gaussian_h_double_prime(x2; FWHM=p.FWHM, H=p.H)

	du[1] = vx1
	du[2] = vz1
	du[3] = - x1f₁ * (p.g + 2 * x1f₂ * vx1^2) / (2 * (1 + x1f₁^2))
	du[4] = (x1f₂ * vx1^2 - 0.5 * p.g * x1f₁^2) / (1 + x1f₁^2)

	du[5] = vx2
	du[6] = vz2
	du[7] = - (0.5 * p.g * x1f₁ * (vx1/vx2) + x2f₁ * x2f₂ * vx2^2) / (1 + x2f₁^2)
	du[8] = (x2f₂ * vx2^2 - 0.5 * p.g * x1f₁ * x2f₁ * (vx1/vx2)) / (1 + x2f₁^2)
end;

# ╔═╡ 19444135-bfcd-45b6-86e3-10a6a9ca36e0
function plot_normal_double_taut_wagon(sol, params; size_w=1000, size_h=1000)
	normal1x = - gaussian_h_double_prime(sol[1, :]; FWHM=params.FWHM, H=params.H) .* gaussian_h_double_prime(sol[1, :]; FWHM=params.FWHM, H=params.H) .* sol[3, :].^2
	normal1x ./= (1 .+ gaussian_h_prime(sol[1, :]; FWHM=params.FWHM, H=params.H).^2)
	
	normal1z = gaussian_h_double_prime(sol[1, :]; FWHM=params.FWHM, H=params.H) .* sol[3, :].^2
	normal1z ./= (1 .+ gaussian_h_prime(sol[1, :]; FWHM=params.FWHM, H=params.H).^2)
	normal1z .+= params.g

	normal1s = sqrt.(normal1x.^2 .+ normal1z.^2)
	
	normal2x = - gaussian_h_double_prime(sol[5, :]; FWHM=params.FWHM, H=params.H) .* gaussian_h_double_prime(sol[5, :]; FWHM=params.FWHM, H=params.H) .* sol[7, :].^2
	normal2x ./= (1 .+ gaussian_h_prime(sol[5, :]; FWHM=params.FWHM, H=params.H).^2)
	
	normal2z = gaussian_h_double_prime(sol[5, :]; FWHM=params.FWHM, H=params.H) .* sol[7, :].^2
	normal2z ./= (1 .+ gaussian_h_prime(sol[5, :]; FWHM=params.FWHM, H=params.H).^2)
	normal2z .+= params.g

	normal2s = sqrt.(normal2x.^2 .+ normal2z.^2)
	
	theme(:juno)
	p1 = plot([0, 0], [-40, 80.0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(p1, [-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(p1, sol[1, :], normal1s; lw=2, label=L"|\mathbf{F_N}|", color="#FE4365")
	plot!(p1, sol[1, :], normal1x; linestyle=:dash, lw=2, label=L"\mathbf{F_{N, x}}", color="#eca25c")
	plot!(p1, sol[1, :], normal1z; linestyle=:dash, lw=2, label=L"\mathbf{F_{N, z}}", color="#3f9778")
	plot!(p1, sol[1, :], sol[2, :]; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(p1, showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!(p1, "Normal Force - Front Wagon", titlelocation=:center)
	xlims!(p1, -50.0, 50.0)
	ylims!(p1, -40, 80)

	p2 = plot([0, 0], [-40, 80.0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(p2, [-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(p2, sol[5, :], normal2s; lw=2, label=L"|\mathbf{F_N}|", color="#FE4365")
	plot!(p2, sol[5, :], normal2x; linestyle=:dash, lw=2, label=L"\mathbf{F_{N, x}}", color="#eca25c")
	plot!(p2, sol[5, :], normal2z; linestyle=:dash, lw=2, label=L"\mathbf{F_{N, z}}", color="#3f9778")
	plot!(p2, sol[5, :], sol[6, :]; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(p2, showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!(p2, "Normal Force - Rear Wagon", titlelocation=:center)
	xlabel!(p2, "x")
	xlims!(p2, -50.0, 50.0)
	ylims!(p2, -40, 80)

	plot(p1, p2; layout = (2, 1))
end;

# ╔═╡ 9400d827-4189-4bc1-b409-1f678c397157
function double_spring_wagon_update!(du, u, p, t)
	x1, z1, vx1, vz1, x2, z2, vx2, vz2 = u

	Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
	S21 = (Delta21 - p.sep) / Delta21

	x1f₁ = gaussian_h_prime(x1; FWHM=p.FWHM, H=p.H)
	x1f₂ = gaussian_h_double_prime(x1; FWHM=p.FWHM, H=p.H)

	x2f₁ = gaussian_h_prime(x2; FWHM=p.FWHM, H=p.H)
	x2f₂ = gaussian_h_double_prime(x2; FWHM=p.FWHM, H=p.H)

	du[1] = vx1
	du[2] = vz1
	
	ax1 = (p.k * S21 * (x2 - x1) - x1f₁ * (x1f₂ * vx1^2 + p.g - p.k * S21 * (z2 - z1))) / (1 + x1f₁^2)
	du[3] = ax1
	du[4] = x1f₂ * vx1^2 + x1f₁ * ax1

	du[5] = vx2
	du[6] = vz2
	
	ax2 = (- p.k * S21 * (x2 - x1) - x2f₁ * (x2f₂ * vx2^2 + p.g + p.k * S21 * (z2 - z1))) / (1 + x2f₁^2) 
	du[7] = ax2
	du[8] = x2f₂ * vx2^2 + x2f₁ * ax2
end;

# ╔═╡ 35e5e538-3f9d-49e1-808d-8d353e02e21c
function n_spring_wagon_update!(du, u, p, t)
	n = p.n
	s0 = p.sep

	for i in 1:n
		idx = 4 * (i - 1)
		
		f₁ = gaussian_h_prime(u[idx + 1]; FWHM=p.FWHM, H=p.H)
		f₂ = gaussian_h_double_prime(u[idx + 1]; FWHM=p.FWHM, H=p.H)

		du[idx + 1] = u[idx + 3]
		du[idx + 2] = u[idx + 4]

		if i == 1
			x1 = u[idx + 1]
			z1 = u[idx + 2]
			vx1 = u[idx + 3]
			vz1 = u[idx + 4]
			
			x2 = u[idx + 5]
			z2 = u[idx + 6]

			Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
			S21 = (Delta21 - s0) / Delta21
			
			ax = (p.k * S21 * (x2 - x1) - f₁ * (f₂ * vx1^2 + p.g - p.k * S21 * (z2 - z1))) / (1 + f₁^2)
			
			du[idx + 3] = ax
			du[idx + 4] = f₂ * vx1^2 + f₁ * ax
		
		elseif i == n
			x1 = u[idx - 3]
			z1 = u[idx - 2]
			
			x2 = u[idx + 1]
			z2 = u[idx + 2]
			vx2 = u[idx + 3]
			vz2 = u[idx + 4]

			Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
			S21 = (Delta21 - s0) / Delta21
			
			ax = (- p.k * S21 * (x2 - x1) - f₁ * (f₂ * vx2^2 + p.g + p.k * S21 * (z2 - z1))) / (1 + f₁^2)
			
			du[idx + 3] = ax
			du[idx + 4] = f₂ * vx2^2 + f₁ * ax
		
		else
			x1 = u[idx - 3]
			z1 = u[idx - 2]
			
			x2 = u[idx + 1]
			z2 = u[idx + 2]
			vx2 = u[idx + 3]
			vz2 = u[idx + 4]

			x3 = u[idx + 5]
			z3 = u[idx + 6]

			Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
			S21 = (Delta21 - s0) / Delta21

			Delta32 = sqrt((x3 - x2)^2 + (z3 - z2)^2)
			S32 = (Delta32 - s0) / Delta32
			
			ax = (- p.k * S21 * (x2 - x1) + p.k * S32 * (x3 - x2) - f₁ * (f₂ * vx2^2 + p.g + p.k * S21 * (z2 - z1) - p.k * S32 * (z3 - z2))) / (1 + f₁^2)
			
			du[idx + 3] = ax
			du[idx + 4] = f₂ * vx2^2 + f₁ * ax
		end
	end
end;

# ╔═╡ 48cde34c-dfa4-4eb6-9b2d-901a2363bacb
function n_spring_wagon_fixed_update!(du, u, p, t)
	n = p.n
	s0 = p.sep / (n - 1)

	for i in 1:n
		idx = 4 * (i - 1)
		
		f₁ = gaussian_h_prime(u[idx + 1]; FWHM=p.FWHM, H=p.H)
		f₂ = gaussian_h_double_prime(u[idx + 1]; FWHM=p.FWHM, H=p.H)

		du[idx + 1] = u[idx + 3]
		du[idx + 2] = u[idx + 4]

		if i == 1
			x1 = u[idx + 1]
			z1 = u[idx + 2]
			vx1 = u[idx + 3]
			vz1 = u[idx + 4]
			
			x2 = u[idx + 5]
			z2 = u[idx + 6]

			Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
			S21 = (Delta21 - s0) / Delta21
			
			ax = (p.k * S21 * (x2 - x1) - f₁ * (f₂ * vx1^2 + p.g - p.k * S21 * (z2 - z1))) / (1 + f₁^2)
			
			du[idx + 3] = ax
			du[idx + 4] = f₂ * vx1^2 + f₁ * ax
		
		elseif i == n
			x1 = u[idx - 3]
			z1 = u[idx - 2]
			
			x2 = u[idx + 1]
			z2 = u[idx + 2]
			vx2 = u[idx + 3]
			vz2 = u[idx + 4]

			Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
			S21 = (Delta21 - s0) / Delta21
			
			ax = (- p.k * S21 * (x2 - x1) - f₁ * (f₂ * vx2^2 + p.g + p.k * S21 * (z2 - z1))) / (1 + f₁^2)
			
			du[idx + 3] = ax
			du[idx + 4] = f₂ * vx2^2 + f₁ * ax
		
		else
			x1 = u[idx - 3]
			z1 = u[idx - 2]
			
			x2 = u[idx + 1]
			z2 = u[idx + 2]
			vx2 = u[idx + 3]
			vz2 = u[idx + 4]

			x3 = u[idx + 5]
			z3 = u[idx + 6]

			Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
			S21 = (Delta21 - s0) / Delta21

			Delta32 = sqrt((x3 - x2)^2 + (z3 - z2)^2)
			S32 = (Delta32 - s0) / Delta32
			
			ax = (- p.k * S21 * (x2 - x1) + p.k * S32 * (x3 - x2) - f₁ * (f₂ * vx2^2 + p.g + p.k * S21 * (z2 - z1) - p.k * S32 * (z3 - z2))) / (1 + f₁^2)
			
			du[idx + 3] = ax
			du[idx + 4] = f₂ * vx2^2 + f₁ * ax
		end
	end
end;

# ╔═╡ a3d63f98-6f4f-4792-b6f6-1299ca1a5a02
function continuum_update(u, p, dt=0.001; tspan=(0.0, 8.0))
	t0, tf = tspan
	steps = round(Int, tf / dt)
	
	n = p.n
	s0 = p.sep / (n - 1)

	a = 0.25 * (s0 / dt)^2

	ts = t0:dt:tf
	nt = length(ts)
	X = Matrix{Float64}(undef, n, nt)

	xprev = u .- p.v0 .* dt
	xcur = u
	xnext = similar(xcur)
	X[:, 1] = xcur

	for j in 1:nt-1
		for i in 1:n
			x = xcur[i]
			f₁ = gaussian_h_prime(x; FWHM=p.FWHM, H=p.H)
			f₂ = gaussian_h_double_prime(x; FWHM=p.FWHM, H=p.H)

			den = 1 + f₁^2

			A = f₁ * f₂ / den
			B = p.g * f₁ / den

			if i == 1
				dx1 = s0 / sqrt(den)
				lap = 2 * (xcur[2] - x) - 2 * dx1
			elseif i == n
				dx1 = s0 / sqrt(den)
				lap = 2 * (xcur[n-1] - x) - 2 * x
			else
				dx1 = (xcur[i+1] - xcur[i-1]) / 2
				lap = xcur[i+1] + xcur[i-1] - 2 * x
			end

			xnext[i] = 2 * x - xprev[i]
					- A * (x - xprev[i])^2 - B * dt^2
					+ a * A * (dt / s0)^2 * dx1^2
					+ a * (dt / s0)^2 * lap
		end
		xprev, xcur, xnext = xcur, xnext, xprev
		X[:, j+1] .= xcur
	end

	return collect(ts), X;
end;

# ╔═╡ 8f84a0e2-0faa-4e24-abd0-ff82fb910557
function plot_normal_continuum(sol_t, sol, params; size_w=1000, size_h=700)
	n = params.n
	s0 = params.sep / (n - 1)
	dt = sol_t[2] - sol_t[1]
	nt = length(sol_t)

	x1 = sol[1, :]
	z1 = gaussian_h(sol[1, :]; FWHM=params.FWHM, H=params.H)

	xdot = similar(sol)
	xdot[:, 2:end-1] .= (sol[:, 3:end] .- sol[:, 1:end-2]) ./ (2 * dt)
	xdot[:, 1] .= (sol[:, 2] .- sol[:, 1]) ./ dt
	xdot[:, end] .= (sol[:, end] .- sol[:, end-1]) ./dt

	FN(i) = @. (gaussian_h_double_prime(sol[i, :]; FWHM=params.FWHM, H=params.H) * xdot[i, :]^2 + g) / (sqrt(1 + gaussian_h_prime(sol[i, :]; FWHM=params.FWHM, H=params.H)^2))
	
	theme(:juno)
	plot([0, 0], [-100, 100]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol[1, :], FN(1); lw=2, label="Normal - Front Wagon", color="#FE4365")
	plot!(sol[n ÷ 2, :], FN(n ÷ 2); lw=2, label="Normal - Middle Wagon", color="#eca25c")
	plot!(sol[end, :], FN(n); lw=2, label="Normal - Rear Wagon", color="#3f9778")
	plot!(x1, z1; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	title!("Normal Force", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-100, 100)
end;

# ╔═╡ 01f26215-273a-4bad-b442-578038df4724
function solve_equations(u0, f, params; tf=20.0)
	tspan = (0.0, 20.0)

	condition(u, t, integrator) = u[1] - 4.0 * integrator.p.FWHM / (2 * sqrt(log(2)))
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)

	prob = ODEProblem(f, u0, tspan, params)
	sol = solve(prob, Tsit5(), callback=cb, dtmax=0.001)

	return sol
end;

# ╔═╡ ee748dde-2466-4b46-816f-938abb217f09
function solve_equations_trisec(u0, f, params; tf=20.0)
	tspan = (0.0, 20.0)

	condition(u, t, integrator) = u[1] - 15.0
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)

	prob = ODEProblem(f, u0, tspan, params)
	sol = solve(prob, Tsit5(), callback=cb, dtmax=0.001)

	return sol
end;

# ╔═╡ e22cd35d-d63f-4740-89b0-2e72b267af9e
function plot_gaussian(FWHM, H; size_w=1000, size_h=700, label=nothing, show_fwhm=false)
	alpha = 2.0 * sqrt(log(2.0)) / FWHM
	xs = range(-50.0, 50.0; length=1000)
	gaussian(x) = H .* exp.(-alpha^2 .* x.^2)

	theme(:juno)
	plot([0, 0], [-10, 60]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(xs, gaussian(xs); size=(size_w, size_h), lw=3, label=label, color="#FE4365")

	if show_fwhm
		plot!([-FWHM/2, FWHM/2], [H/2, H/2]; label=nothing, linestyle=:dash, color="#f2f3a2")
		plot!([-FWHM/2, -FWHM/2], [0, H/2]; label=nothing, linestyle=:dash, color="#f2f3a2")
		plot!([FWHM/2, FWHM/2], [0, H/2]; label=nothing, linestyle=:dash, color="#f2f3a2")
	end
	
	plot!(showaxis=:no,ticks=nothing)
	xlims!(-50.0, 50.0)
	ylims!(-10.0, 60.0)
end;

# ╔═╡ d8e73ce1-4444-4270-87be-5e5a619b49b1
function plot_trisectrix(b; size_w=1000, size_h=700, label=nothing)
	u = collect(-100:0.05:100)
	xs = b .* u .* (u.^2 .- 3) ./ (u.^2 .+ 1)
	zs = - b .* (u.^2 .- 3) ./ (u.^2 .+ 1)
	# println(b)
	
	theme(:juno)
	plot([0, 0], [-10, 60]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(xs, zs; size=(size_w, size_h), lw=3, label=label, color="#FE4365")
	
	plot!(showaxis=:no,ticks=nothing)
	xlims!(-50.0, 50.0)
	ylims!(-10.0, 60.0)
end;

# ╔═╡ a8c6c34c-b6d6-4f77-b6cd-ebe9f8776271
function calculate_spring_normal(sol, params, n, sep_val=1.0; sep=false, n_wagon=false, fixed=false)
	if sep
		s0 = params.sep
	else
		s0 = sep_val
	end

	if n_wagon
		n0 = params.n
	else
		n0 = 2
	end

	if fixed
		s0 /= (n0 - 1)
	end
	
	normal = similar(sol[1, :])
	idx = 4 * (n - 1)

	f₁ = gaussian_h_prime(sol[idx + 1, :]; FWHM=params.FWHM, H=params.H)
	f₂ = gaussian_h_double_prime(sol[idx + 1, :]; FWHM=params.FWHM, H=params.H)
	
	if n == 1
		x1 = sol[idx + 1, :]
		z1 = sol[idx + 2, :]
		vx1 = sol[idx + 3, :]
		vz1 = sol[idx + 4, :]
		
		x2 = sol[idx + 5, :]
		z2 = sol[idx + 6, :]

		Delta21 = similar(x2)
		S21 = similar(x2)
		
		@. Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
		@. S21 = (Delta21 - s0) / Delta21

		@. normal = (params.g + f₂ * vx1^2 + params.k * S21 * (f₁ * (x2 - x1) - (z2 - z1))) / (1 + f₁^2)
	
	elseif n == n0
		x1 = sol[idx - 3, :]
		z1 = sol[idx - 2, :]
		
		x2 = sol[idx + 1, :]
		z2 = sol[idx + 2, :]
		vx2 = sol[idx + 3, :]
		vz2 = sol[idx + 4, :]

		Delta21 = similar(x2)
		S21 = similar(x2)

		@. Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
		@. S21 = (Delta21 - s0) / Delta21
		
		@. normal = (params.g + f₂ * vx2^2 - params.k * S21 * (f₁ * (x2 - x1) - (z2 - z1))) / (1 + f₁^2)
	
	else
		x1 = sol[idx - 3, :]
		z1 = sol[idx - 2, :]
		
		x2 = sol[idx + 1, :]
		z2 = sol[idx + 2, :]
		vx2 = sol[idx + 3, :]
		vz2 = sol[idx + 4, :]

		x3 = sol[idx + 5, :]
		z3 = sol[idx + 6, :]
		
		Delta21 = similar(x2)
		S21 = similar(x2)

		Delta32 = similar(x2)
		S32 = similar(x2)

		@. Delta21 = sqrt((x2 - x1)^2 + (z2 - z1)^2)
		@. S21 = (Delta21 - s0) / Delta21

		@. Delta32 = sqrt((x3 - x2)^2 + (z3 - z2)^2)
		@. S32 = (Delta32 - s0) / Delta32
		
		@. normal = (params.g + f₂ * vx2^2 - params.k * S21 * (f₁ * (x2 - x1) - (z2 - z1)) + params.k * S32 * (f₁ * (x3 - x2) - (z3 - z2))) / (1 + f₁^2)
	end

	return normal
end;

# ╔═╡ ff92cc16-fcda-4032-bd5f-52007b50b204
function plot_normal_double_spring_wagon(sol, params; size_w=1000, size_h=700)
	normal1s = calculate_spring_normal(sol, params, 1)
	normal2s = calculate_spring_normal(sol, params, 2)
	
	theme(:juno)
	plot([0, 0], [-40, 80]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol[1, :], normal1s; lw=2, label="Normal - Front Wagon", color="#FE4365")
	plot!(sol[5, :], normal2s; lw=2, label="Normal - Rear Wagon", color="#eca25c")
	plot!(sol[1, :], sol[2, :]; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legend=:bottomright, legendfontsize=12)
	title!("Normal Force", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-40, 80)
end;

# ╔═╡ aa71d89e-2190-4b39-bf1d-c7a680c3d1d9
function plot_normal_n_spring_wagon(sol, params; size_w=1000, size_h=700, fixed=false)
	n = params.n
	middle = floor(Int, n/2)

	normal1s = calculate_spring_normal(sol, params, 1; sep=true, n_wagon=true, fixed=fixed)
	normal2s = calculate_spring_normal(sol, params, middle; sep=true, n_wagon=true, fixed=fixed)
	normal3s = calculate_spring_normal(sol, params, n; sep=true, n_wagon=true, fixed=fixed)

	theme(:juno)
	plot([0, 0], [-40, 80]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!([-50, 50], [0, 0]; linestyle=:dot, lw=2, label=nothing, color="#a1a1a1")
	plot!(sol[1, :], normal1s; lw=2, label="Normal - Front Wagon", color="#FE4365")
	plot!(sol[4 * middle + 1, :], normal2s; lw=2, label="Normal - Middle Wagon", color="#eca25c")
	plot!(sol[end - 3, :], normal3s; lw=2, label="Normal - Rear Wagon", color="#3f9778")
	plot!(sol[1, :], sol[2, :]; lw=1, linestyle=:dash, label=nothing, color="#f2f3a2")
	plot!(showaxis=:no, ticks=nothing, size=(size_w, size_h), legendfontsize=12)
	title!("Normal Force", titlelocation=:center)
	xlabel!("x")
	xlims!(-50.0, 50.0)
	ylims!(-40, 80)
end;

# ╔═╡ 90c88c7d-ce7d-4db8-baeb-2cd3698b239a
md"""
## Sliders
"""

# ╔═╡ e62a5f89-47eb-4841-bde8-2b018eae0fba
function gaussian_sliders(;show_val=false)
	return PlutoUI.combine() do Child
		md"""
		#### Gaussian parameters
		``
		z = h(x) = H \exp\left(-\alpha^2 x^2\right) \quad ; \quad \alpha = \frac{2\sqrt{\ln 2}}{\text{FWHM}}\qquad``
		FWHM: $(Child("FWHM", Slider(1.0:1.0:50.0; default=20.0, show_value=show_val)))
		Height: $(Child("H", Slider(1.0:1.0:50.0; default=30.0, show_value=show_val)))
		"""
	end
end;

# ╔═╡ 910d4b85-999e-4c7c-87b2-7cb3f294dbea
@bind blind_params gaussian_sliders()

# ╔═╡ 41318e0f-f449-4605-8def-bd0ec59c4f39
PlutoCSSMacros.css_horizontal(track_gaussian_html, plot_gaussian(blind_params.FWHM, blind_params.H; show_fwhm=show_fwhm); width=[45, 55])

# ╔═╡ 94ee28f9-c844-4fef-a725-93741318bed0
begin
	md"""
	$(@bind single_wagon_params gaussian_sliders())
	---
	"""
end

# ╔═╡ d24d621e-2350-429b-8def-dade2b04517f
function trisectrix_sliders(;show_val=false)
	return PlutoUI.combine() do Child
		md"""
		#### Trisectrix parameters
		``x(u) = \dfrac{b u (u^2 - 3)}{u^2 + 1} ; \qquad \qquad z(u) = - \dfrac{b (u^2 - 3)}{u^2 + 1}``
		
		b: $(Child("b", Slider(1.0:0.1:20.0; default=5.0, show_value=show_val)))
		"""
	end
end;

# ╔═╡ 6f8ce07d-7128-498a-a289-170b64e446af
@bind trisec_params trisectrix_sliders(show_val=true)

# ╔═╡ edcb0089-33b9-44aa-bea6-b633840dbe66
PlutoCSSMacros.css_horizontal(track_trisectrix_html, plot_trisectrix(trisec_params.b); width=[45, 55])

# ╔═╡ dbec0cc5-61a2-4dec-bdce-663862ceb90f
function x_slider(; show_val=false)
	return PlutoUI.combine() do Child
		md"""
		``x_0``: $(Child("x", Slider(-50.0:1.0:50.0; default=-40.0, show_value=show_val)))
		"""
	end
end;

# ╔═╡ 169d7dcc-ce3e-4483-b563-af302dd82dd4
md"""
### Position of the wagon:

$(@bind x_single_wagon x_slider())

show forces: $(@bind show_forces Switch()) ``\qquad`` show tangent: $(@bind show_tangent Switch())
"""

# ╔═╡ ef49182e-eda4-4217-9a88-a3d5f4a1eade
function plot_forces(x_single_wagon, single_wagon_params; show_forces=false, show_tangent=show_tangent, size_w=1000, size_h=700)
	FWHM = single_wagon_params.FWHM
	H = single_wagon_params.H
	x0 = x_single_wagon.x
	z0 = gaussian_h(x0, FWHM=FWHM, H=H)

	mult = 10.0
	
	plt = plot_gaussian(FWHM, H; size_w=size_w, size_h=size_h)
	tangent_angle = atan(gaussian_h_prime(x0; FWHM=FWHM, H=H))
	
	if show_forces
		normal_x = [x0, x0 - mult * sin(tangent_angle)]
		normal_z = [z0, z0 + mult * cos(tangent_angle)]
		
		plot!(plt, normal_x, normal_z; arrow=:closed, lw=2, color="darkcyan", label="Normal")
		plot!(plt, [x0, x0], [z0, z0-mult]; arrow=:closed, lw=2, color="goldenrod", label="Gravity")
	end

	if show_tangent
		tangent_x = [x0, x0 + mult * cos(tangent_angle)]
		tangent_z = [z0, z0 + mult * sin(tangent_angle)]
		
		plot!(plt, tangent_x, tangent_z; arrow=:closed, lw=2, color="skyblue", label="Direction of velocity")
	end
	
	scatter!(plt, [x0], [z0]; markersize=10, color="cyan", label="wagon")
	plot!(plt, legend=:topleft)
end;

# ╔═╡ dfa26e29-10d5-46c4-8f5f-dc19e8133fd3
@PlutoCSSMacros.css_center plot_forces(x_single_wagon, single_wagon_params; show_forces=show_forces, show_tangent=show_tangent, size_w=600, size_h=420)

# ╔═╡ a1f99534-45b0-4ebf-bc00-374de1c0a8be
function plot_forces_double(x_double_wagon, params; show_forces=false, show_tangent=show_tangent, size_w=1000, size_h=700)
	FWHM = params.FWHM
	H = params.H
	x1 = x_double_wagon.x
	z1 = gaussian_h(x1; FWHM=FWHM, H=H)

	sep = 10.0

	fx(x2) = (x2 - x1)^2 + (gaussian_h(x2; FWHM=FWHM, H=H) - gaussian_h(x1; FWHM=FWHM, H=H))^2 - sep^2
	x2 = find_zero(fx, (x1-sep-1, x1))
	
	z2 = gaussian_h(x2; FWHM = FWHM, H=H)
	
	mult = 10.0
	
	plt = plot_gaussian(FWHM, H; size_w=size_w, size_h=size_h)
	tangent_angle1 = atan(gaussian_h_prime(x1; FWHM=FWHM, H=H))
	tangent_angle2 = atan(gaussian_h_prime(x2; FWHM=FWHM, H=H))

	if show_forces
		normal_x1 = [x1, x1 - mult * sin(tangent_angle1)]
		normal_z1 = [z1, z1 + mult * cos(tangent_angle1)]

		normal_x2 = [x2, x2 - mult * sin(tangent_angle2)]
		normal_z2 = [z2, z2 + mult * cos(tangent_angle2)]
		
		plot!(plt, normal_x1, normal_z1; arrow=:closed, lw=2, color="darkcyan", label="Normal")
		plot!(plt, normal_x2, normal_z2; arrow=:closed, lw=2, color="darkcyan", label=nothing)
		plot!(plt, [x1, x1], [z1, z1-mult]; arrow=:closed, lw=2, color="goldenrod", label="Gravity")
		plot!(plt, [x2, x2], [z2, z2-mult]; arrow=:closed, lw=2, color="goldenrod", label=nothing)
	end

	if show_tangent
		tangent_x1 = [x1, x1 + mult * cos(tangent_angle1)]
		tangent_z1 = [z1, z1 + mult * sin(tangent_angle1)]

		tangent_x2 = [x2, x2 + mult * cos(tangent_angle2)]
		tangent_z2 = [z2, z2 + mult * sin(tangent_angle2)]
		
		plot!(plt, tangent_x1, tangent_z1; arrow=:closed, lw=2, color="skyblue", label="Direction of velocity")
		plot!(plt, tangent_x2, tangent_z2; arrow=:closed, lw=2, color="skyblue", label=nothing)
	end
	
	scatter!(plt, [x1, x2], [z1, z2]; markersize=10, color="cyan", label="wagon")
	plot!(plt, legend=:topleft)
end;

# ╔═╡ 10bdc6ae-eb28-4a4c-bdbd-3ebe883e1a3c
md"""
### Position of the front wagon:

$(@bind x_double_taut_wagon x_slider())

show forces: $(@bind show_forces_double Switch()) ``\qquad`` show tangent: $(@bind show_tangent_double Switch())
"""

# ╔═╡ 001027a3-3390-412e-a168-40565c203809
@PlutoCSSMacros.css_center plot_forces_double(x_double_taut_wagon, params; show_forces=show_forces_double, show_tangent=show_tangent_double, size_w=600, size_h=420)

# ╔═╡ 38e6b553-63be-462a-88c4-eea9854ba73e
function vx_slider(; show_val=false, v_min=0.8, v_max=1.2, default=1.018)
	return PlutoUI.combine() do Child
		md"""
		``v_{x,0}``: $(Child("vx", Scrubbable(v_min:0.001:v_max; format=".3f", default=default))) ``v_{min}``
		"""
	end
end;

# ╔═╡ c7fc0b70-9203-46c7-83f0-3f331aea129f
md"""
### Initial conditions of the wagon:
``m = 1; \qquad \qquad x_0 = -\dfrac{4}{\alpha}; \qquad \qquad v_{\text{min}} = \sqrt{2 \, g \, H}; \qquad \qquad N = - \dfrac{\ddot{x}}{\sin \theta} = \dfrac{g + h'' \, \dot{x}^2}{\sqrt{1 + h'^2}}``

$(@bind vx_single_wagon vx_slider())
"""

# ╔═╡ 4f62c0af-9077-4ecd-8d88-422c3f13b832
begin
	vx_single = vx_single_wagon.vx * sqrt(2 * params.g * params.H)
	u0_single_wagon = [
		x0,
		gaussian_h(x0; FWHM=params.FWHM, H=params.H),
		vx_single,
		vx_single * gaussian_h_prime(x0; FWHM=params.FWHM, H=params.H)
	]

	sol_single_wagon = solve_equations(u0_single_wagon, single_wagon_update!, params)
end;

# ╔═╡ 16d57df6-161e-446b-a329-09c67fbe5b20
PlutoCSSMacros.css_horizontal(plot_energy_single_wagon(sol_single_wagon, params; size_w=600, size_h=450), plot_normal_single_wagon(sol_single_wagon, params; size_w=600, size_h=450))

# ╔═╡ 385bfa88-e4da-4549-9c6b-eb7ee5008ff5
md"""
### Initial conditions of the wagons:
``x_{1, 0} = -\dfrac{4}{\alpha}; \qquad \qquad x_{2, 0} = x_{1, 0} - s_0; \qquad \qquad v_{1, 0} = v_{2, 0} = v_{x, 0}``

$(@bind vx_double_free_wagon vx_slider())
"""

# ╔═╡ aebdbc62-c568-4f30-af3c-2d287d477177
begin
	vx_double_free = vx_double_free_wagon.vx * sqrt(2 * params.g * params.H)
	
	u0_double_free_wagon = [
		x0,
		gaussian_h(x0; FWHM=params.FWHM, H=params.H),
		vx_double_free,
		vx_double_free * gaussian_h_prime(x0; FWHM=params.FWHM, H=params.H),
		x0 - 1,
		gaussian_h(x0 - 1; FWHM=params.FWHM, H=params.H),
		vx_double_free,
		vx_double_free * gaussian_h_prime(x0 - 1; FWHM=params.FWHM, H=params.H),
	]

	sol_double_free_wagon = solve_equations(u0_double_free_wagon, double_free_wagon_update!, params)
end;

# ╔═╡ 01ae5c23-5ee7-404a-bf00-8b7d126bb257
PlutoCSSMacros.css_horizontal(plot_separation_double_free_wagon(sol_double_free_wagon, params; size_w=600, size_h=450), plot_normal_double_free_wagon(sol_double_free_wagon, params; size_w=600, size_h=450))

# ╔═╡ c0518b72-3df1-402f-a7e5-ebb54332ad5e
md"""
### Initial conditions of the wagons:
``x_{1, 0} = -\dfrac{4}{\alpha}; \qquad \qquad x_{2, 0} = x_{1, 0} - s_0; \qquad \qquad v_{1, 0} = v_{x, 0}``

$(@bind vx_double_taut_wagon vx_slider(; v_min=0.7, v_max=1.2, default=0.8))
``s_0:`` $(@bind separation Slider(0.5:0.1:2.5; default=1.0))
"""

# ╔═╡ 85ba38eb-33c4-4596-a5fa-3691b77d4d7e
begin
	vx_double_taut = vx_double_taut_wagon.vx * sqrt(2 * params.g * params.H)
	
	u0_double_taut_wagon = [
		x0,
		gaussian_h(x0; FWHM=params.FWHM, H=params.H),
		vx_double_taut,
		vx_double_taut * gaussian_h_prime(x0; FWHM=params.FWHM, H=params.H),
		
		x0 - separation,
		gaussian_h(x0 - separation; FWHM=params.FWHM, H=params.H),
		vx_double_taut * sqrt((1 + gaussian_h_prime(x0; FWHM=params.FWHM, H=params.H)^2) / (1 + gaussian_h_prime(x0 - separation; FWHM=params.FWHM, H=params.H)^2)),
		vx_double_taut * sqrt((1 + gaussian_h_prime(x0; FWHM=params.FWHM, H=params.H)^2) / (1 + gaussian_h_prime(x0 - separation; FWHM=params.FWHM, H=params.H)^2)) * gaussian_h_prime(x0 - separation; FWHM=params.FWHM, H=params.H)
	]

	sol_double_taut_wagon = solve_equations(u0_double_taut_wagon, double_taut_wagon_update!, params)
end;

# ╔═╡ bb729bad-2aa6-4080-80cf-24469010a3fc
PlutoCSSMacros.css_horizontal(plot_separation_double_taut_wagon(sol_double_taut_wagon, params; size_w=600, size_h=450), plot_normal_double_taut_wagon(sol_double_taut_wagon, params; size_w=600, size_h=600))

# ╔═╡ ece97e65-5310-4e4a-aa2d-0aac3a862455
md"""
### Initial conditions of the wagons:
``x_{1, 0} = -\dfrac{4}{\alpha}; \qquad \qquad x_{2, 0} = x_{1, 0} - s_0; \qquad \qquad v_{1, 0}= v_{x, 0}``

``\frac{k}{m}:`` $(@bind hardness Slider(1:1:30; default=20, show_value=true)) ``\frac{g}{s_0}``

$(@bind vx_double_spring_wagon vx_slider(;default=1.080))
"""

# ╔═╡ 57f88144-537f-4c90-86f1-ee0481e00e4c
begin
	vx_double_spring = vx_double_spring_wagon.vx * sqrt(2 * params.g * params.H)
	
	u0_double_spring_wagon = [
		x0,
		gaussian_h(x0; FWHM=params.FWHM, H=params.H),
		vx_double_spring,
		vx_double_spring * gaussian_h_prime(x0; FWHM=params.FWHM, H=params.H),
		
		x0 - 1,
		gaussian_h(x0 - 1; FWHM=params.FWHM, H=params.H),
		vx_double_spring,
		vx_double_spring * gaussian_h_prime(x0 - 1; FWHM=params.FWHM, H=params.H)
	]

	double_spring_params = (
		FWHM = 20.0,
		H = 30.0,
		g = 10.0,
		sep = 1.0,
		k = hardness * 10.0 / 1.0 
	)

	sol_double_spring_wagon = solve_equations(u0_double_spring_wagon, double_spring_wagon_update!, double_spring_params)
end;

# ╔═╡ 107bc604-4eac-4755-852e-7ec0420a2a8a
PlutoCSSMacros.css_horizontal(plot_normal_double_spring_wagon(sol_double_spring_wagon, double_spring_params; size_w=600, size_h=450), plot_velocity_double_spring_wagon(sol_double_spring_wagon, double_spring_params; size_w=600, size_h=450, time=false))

# ╔═╡ 971eaae6-057c-46cf-9258-5168dcc7b0a7
PlutoCSSMacros.css_horizontal(plot_separation_double_spring_wagon(sol_double_spring_wagon, double_spring_params; size_w=600, size_h=450), plot_velocity_double_spring_wagon(sol_double_spring_wagon, double_spring_params; size_w=600, size_h=450))

# ╔═╡ 5fb5b59d-97a9-4292-95b1-b957d39823af
md"""
### Initial conditions of the wagons:
``x_{1, 0} = -\dfrac{4}{\alpha}; \qquad \qquad x_{2, 0} = x_{1, 0} - s_0; \qquad \qquad v_{1, 0}= v_{x, 0}``

``\frac{k}{m}:`` $(@bind hardness_n_wagon Slider(10:1:30; default=20, show_value=true)) ``\frac{g}{s_0}`` ``\qquad \qquad N:`` $(@bind n_wagon Slider(3:1:15; default=3, show_value=true))

$(@bind vx_n_spring_wagon vx_slider(;default=1.018))
"""

# ╔═╡ 1c133ded-ff67-487a-85a9-ec24084f0dc4
begin
	vx_n_spring = vx_n_spring_wagon.vx * sqrt(2 * params.g * params.H)

	u0_n_spring_wagon = Float64[]
	for i in 1:n_wagon
		push!(u0_n_spring_wagon,
			  x0 - (i - 1) * 1.0,
			  gaussian_h(x0 - (i - 1) * 1.0; FWHM=params.FWHM, H=params.H),
			  vx_n_spring,
			  vx_n_spring * gaussian_h_prime(x0 - (i - 1) * 1.0; FWHM=params.FWHM, H=params.H)
			 )
	end

	n_spring_params = (
		FWHM = 20.0,
		H = 30.0,
		g = 10.0,
		sep = 1.0,
		k = hardness_n_wagon * 10.0 / 1.0,
		n = n_wagon
	)

	sol_n_spring_wagon = solve_equations(u0_n_spring_wagon, n_spring_wagon_update!, n_spring_params)
end;

# ╔═╡ 82e023e7-e03c-4c6a-a34a-aa690d15539f
PlutoCSSMacros.css_horizontal(plot_normal_n_spring_wagon(sol_n_spring_wagon, n_spring_params; size_w=600, size_h=450), plot_velocity_n_spring_wagon(sol_n_spring_wagon, n_spring_params; size_w=600, size_h=450, time=false))

# ╔═╡ 0dacebcc-017f-4bb4-92cb-1749fd574f27
PlutoCSSMacros.css_horizontal(plot_separation_n_spring_wagon(sol_n_spring_wagon, n_spring_params; size_w=600, size_h=450), plot_velocity_n_spring_wagon(sol_n_spring_wagon, n_spring_params; size_w=600, size_h=450))

# ╔═╡ 9f631ba6-bebb-46e8-bd5c-596f5629ff27
md"""
### Initial conditions of the wagons:
``x_{1, 0} = -\dfrac{4}{\alpha}; \qquad \qquad x_{2, 0} = x_{1, 0} - \frac{L}{N - 1}; \qquad \qquad v_{1, 0}= v_{x, 0}``

``\frac{k}{m}:`` $(@bind hardness_n_wagon_fixed Slider(10:1:30; default=20, show_value=true)) ``\frac{g}{s_0}`` ``\qquad \qquad N:`` $(@bind n_wagon_fixed Slider(3:1:25; default=3, show_value=true))

$(@bind vx_n_spring_wagon_fixed vx_slider(;default=1.018))
"""

# ╔═╡ c02b530f-b30d-4abd-8498-8b2689078cb7
begin
	vx_n_spring_fixed = vx_n_spring_wagon_fixed.vx * sqrt(2 * params.g * params.H)

	u0_n_spring_wagon_fixed = Float64[]
	for i in 1:n_wagon_fixed
		push!(u0_n_spring_wagon_fixed,
			  x0 - (i - 1) * 5.0 / (n_wagon_fixed - 1),
			  gaussian_h(x0 - (i - 1) * 5.0 / (n_wagon_fixed - 1); FWHM=params.FWHM, H=params.H),
			  vx_n_spring_fixed,
			  vx_n_spring_fixed * gaussian_h_prime(x0 - (i - 1) * 5.0 / (n_wagon_fixed - 1); FWHM=params.FWHM, H=params.H)
			 )
	end

	n_spring_fixed_params = (
		FWHM = 20.0,
		H = 30.0,
		g = 10.0,
		sep = 5.0,
		k = hardness_n_wagon_fixed * 10.0 * n_wagon_fixed / 5.0,
		n = n_wagon_fixed
	)

	sol_n_spring_wagon_fixed = solve_equations(u0_n_spring_wagon_fixed, n_spring_wagon_fixed_update!, n_spring_fixed_params)
end;

# ╔═╡ 127b4a43-9537-42f5-a30f-10858012ceef
PlutoCSSMacros.css_horizontal(plot_normal_n_spring_wagon(sol_n_spring_wagon_fixed, n_spring_fixed_params; size_w=600, size_h=450, fixed=true), plot_velocity_n_spring_wagon(sol_n_spring_wagon_fixed, n_spring_fixed_params; size_w=600, size_h=450, time=false))

# ╔═╡ dd6bec64-7b56-47a4-ade2-790008285aaf
PlutoCSSMacros.css_horizontal(plot_separation_n_spring_wagon(sol_n_spring_wagon_fixed, n_spring_fixed_params; size_w=600, size_h=450, fixed=true), plot_velocity_n_spring_wagon(sol_n_spring_wagon_fixed, n_spring_fixed_params; size_w=600, size_h=450))

# ╔═╡ a6dacc52-83ac-4f80-96de-ac1521246592
md"""
### Initial conditions of the wagons:
``x_{0} = -\dfrac{4}{\alpha}; \qquad \qquad s_0 = \dfrac{L}{N - 1}, N \to \infty; \qquad \qquad v_0 = v_{x, 0}``

$(@bind vx_continuum_wagon vx_slider(;v_min=0.5, v_max=0.9, default=0.6))
"""

# ╔═╡ eb7eeff2-5314-44d5-b10a-e07b1c68a600
begin
	vx_continuum = vx_continuum_wagon.vx * sqrt(2 * params.g * params.H)
	n_continuum = 51

	u0_continuum = Float64[]
	for i in 1:n_continuum
		push!(u0_continuum,
			  x0 - (i - 1) * 10.0 / (n_continuum - 1)
			 )
	end

	n_continuum_params = (
		FWHM = 20.0,
		H = 30.0,
		g = 10.0,
		sep = 5.0,
		v0 = vx_continuum,
		n = n_continuum
	)

	sol_continuum_t, sol_continuum = continuum_update(u0_continuum, n_continuum_params)
end;

# ╔═╡ 15837e7d-b35d-4e11-bb4c-e24047771e46
# PlutoCSSMacros.css_horizontal(plot_trajectory_continuum(sol_continuum_t, sol_continuum, n_continuum_params; size_w=600, size_h=450, percent_t=percent_trajectory), plot_normal_continuum(sol_continuum_t, sol_continuum, n_continuum_params; size_w=600, size_h=450))

@PlutoCSSMacros.css_center plot_trajectory_continuum(sol_continuum_t, sol_continuum, n_continuum_params; size_w=600, size_h=450, percent_t=percent_trajectory)

# ╔═╡ 0d7c01fa-5e3f-482f-aa4b-1824ef00a275
md"""
### Initial conditions of the wagon:
``m = 1; \qquad \qquad u_0 = -15; \qquad \qquad v_{\text{min}} = \sqrt{2 \, g \, (z_{\text{max}} - z_{\text{min}})}``

``b:`` $(@bind b_trisec Slider(1.0:0.1:20.0; default=5))

$(@bind vx_single_wagon_trisec vx_slider())
"""

# ╔═╡ 816cb29e-b8a6-4579-a755-93709b22e967
begin
	single_wagon_trisec_params = (
		g = 10.0,
		H = 4 * b_trisec,
		b = b_trisec
	)
	
	vx_single_trisec = vx_single_wagon_trisec.vx * sqrt(2 * single_wagon_trisec_params.g * single_wagon_trisec_params.H)
	u0_single_wagon_trisec = [-15.0, vx_single_trisec / trisec_dx(-15.0, single_wagon_trisec_params.b)]
	
	sol_single_wagon_trisec = solve_equations_trisec(u0_single_wagon_trisec, single_wagon_trisectrix_update!, single_wagon_trisec_params)
end;

# ╔═╡ b6f990b7-81ec-45e6-9662-a85da2287f23
PlutoCSSMacros.css_horizontal(plot_energy_single_wagon_trisec(sol_single_wagon_trisec, single_wagon_trisec_params; size_w=600, size_h=450), plot_normal_single_wagon_trisec(sol_single_wagon_trisec, single_wagon_trisec_params; size_w=600, size_h=450))

# ╔═╡ b25185c9-3e1b-4cc9-9b24-54ac7492a0a1
md"""
# Extras
"""

# ╔═╡ Cell order:
# ╟─ff3f2dee-a847-4df5-bf8c-6fc7c22ebe70
# ╟─ab865418-5e3a-4062-ad56-e4df175f0128
# ╟─8cbf05ba-b6a4-4e9c-ac38-b7a3963ca1b5
# ╟─1995dd8a-5ee0-4c26-8dfc-95135f36bf50
# ╟─09a2545f-7aa9-4a9c-a215-e1e710ccd99d
# ╟─55923503-d43f-48fb-b0e9-d44cea7043f2
# ╟─81b8bfb2-5441-4e16-9e04-739c44c72a71
# ╟─9bb46c1f-2b04-445b-a245-e2d12f14f94d
# ╟─0a02ba54-5256-4b0e-9312-cabc6e8cce45
# ╟─d0151128-78e7-43a6-8b5c-e595e30dd25d
# ╟─f738b9b4-ee60-4a2c-9cc2-a15beb7cbcc6
# ╟─71b1e8fc-e5dc-41a1-8220-a36036bfd3ba
# ╟─851c6405-c9de-49a1-9caa-1f2c66423012
# ╟─69bad9fd-c9c5-4699-93d1-103505692e63
# ╟─ba2126d7-43bc-4878-9f4a-fd6dfc2749b6
# ╟─517d9166-36e2-4fb6-b991-404b8023f811
# ╟─fd72801e-3131-4950-8185-bc382d1e4618
# ╟─3a2f3ffc-189e-4450-966b-73173f8df444
# ╟─a835de6d-e0ed-481a-a5d1-db0976e251c0
# ╟─910d4b85-999e-4c7c-87b2-7cb3f294dbea
# ╟─a1c7d6ac-bb2e-480c-8469-44bbb1c99e60
# ╟─41318e0f-f449-4605-8def-bd0ec59c4f39
# ╟─735e1043-20c8-4e46-a576-80eab7c8707e
# ╟─1bf8a803-bb08-4e55-b014-991ec14c6824
# ╟─94ee28f9-c844-4fef-a725-93741318bed0
# ╟─169d7dcc-ce3e-4483-b563-af302dd82dd4
# ╟─1ce1a72c-85de-41b8-9de7-9c3548f94791
# ╟─dfa26e29-10d5-46c4-8f5f-dc19e8133fd3
# ╟─ef49182e-eda4-4217-9a88-a3d5f4a1eade
# ╟─6fa60588-963d-4b6b-95ff-37b5ba75b20c
# ╟─d5ce4638-3890-49ea-91f9-15b54c6e96e1
# ╟─31bb8c30-714c-4b72-a0fc-b9ae167a5fee
# ╟─12fe1903-e002-4749-be65-6554609320c7
# ╟─cef6e64f-d9fd-4842-a7be-ff46f3787a8d
# ╟─63d25314-eee0-408f-b3a1-eefafa81a953
# ╟─51b15b08-75f2-4c22-b515-c0b1cb359a5f
# ╟─ecf6d1bb-e076-4481-a8e1-59afeec5609f
# ╟─af81083b-c29a-419b-90a7-3223a198fe87
# ╟─7a2df5f6-5ab9-4759-96ad-e60e44200fb1
# ╟─a93a5477-2586-4b12-8841-77a1ba3a284c
# ╟─bd3fa162-19b3-4ca8-84dc-4a7e662c0000
# ╟─c72f6a42-4215-4eee-a9c8-97cc62855384
# ╟─313af5c9-ee74-4655-ae4e-240ef9dc6a22
# ╟─7302dd08-8601-451e-91f5-5411dd670750
# ╟─69f15733-5946-40ef-9c37-15b5f4c6b73b
# ╟─c7fc0b70-9203-46c7-83f0-3f331aea129f
# ╟─4f62c0af-9077-4ecd-8d88-422c3f13b832
# ╟─16d57df6-161e-446b-a329-09c67fbe5b20
# ╟─b7b7ff06-32e9-45f4-b043-b32a1c6f3f01
# ╟─ac303828-ea10-4266-a44a-826623740958
# ╟─da9420e4-1ddc-400a-8993-b6fee5eda3ae
# ╟─385bfa88-e4da-4549-9c6b-eb7ee5008ff5
# ╟─aebdbc62-c568-4f30-af3c-2d287d477177
# ╟─01ae5c23-5ee7-404a-bf00-8b7d126bb257
# ╟─ddbf78aa-29ce-41e2-ba61-45aab0d9c942
# ╟─4163a1cf-8e8e-4337-9796-a9d7c6cf8608
# ╟─8b2b545c-5ed0-44be-a881-63096994ba62
# ╟─79eb156c-7fd4-4588-acbc-cbdbeec95a18
# ╟─10bdc6ae-eb28-4a4c-bdbd-3ebe883e1a3c
# ╟─001027a3-3390-412e-a168-40565c203809
# ╟─a1f99534-45b0-4ebf-bc00-374de1c0a8be
# ╟─0e78beae-b936-430c-8f08-11c6a48fb884
# ╟─afc900d4-557c-4d1e-91cd-ae89ff027291
# ╟─524258a3-440e-44f3-85b3-c28a58903836
# ╟─cdccc6e1-c805-4d07-95fa-7b02bd51b7c9
# ╟─39776d7a-1ba8-4821-9211-978e0ade0554
# ╟─439eda1f-fdcb-4cb2-b06c-275ed6f1b16c
# ╟─ad911bfc-3f35-4173-a849-765a0ee2c88d
# ╟─85f5bb91-2b5c-4a36-a711-738264410ae0
# ╟─c0518b72-3df1-402f-a7e5-ebb54332ad5e
# ╟─85ba38eb-33c4-4596-a5fa-3691b77d4d7e
# ╟─bb729bad-2aa6-4080-80cf-24469010a3fc
# ╟─078dbc7f-e3c4-44fa-aa86-4de02574a0a5
# ╟─19444135-bfcd-45b6-86e3-10a6a9ca36e0
# ╟─1d4d6278-8e6b-4d93-a0eb-586202d4e4c1
# ╟─149b825b-3ea7-4aeb-9d9e-94291ddd0bf7
# ╟─9d6624a4-3142-41bb-aebf-1ebbf66a20a7
# ╟─1ce2148a-49a8-4fb5-822d-e4bf4b53fdb7
# ╟─36fd7fe2-e34b-4118-9bb6-a797ec963bfc
# ╟─72dfbb63-87c6-4218-9237-0a1b8e56493a
# ╟─730c2c53-bc09-4783-a40a-2bdc5c21c25c
# ╟─ece97e65-5310-4e4a-aa2d-0aac3a862455
# ╟─57f88144-537f-4c90-86f1-ee0481e00e4c
# ╟─107bc604-4eac-4755-852e-7ec0420a2a8a
# ╟─ff92cc16-fcda-4032-bd5f-52007b50b204
# ╟─9400d827-4189-4bc1-b409-1f678c397157
# ╟─158ad312-bbf0-492b-aede-23ae67ee0198
# ╟─2b82e765-6017-4d92-8c0d-b4aab70a9110
# ╟─971eaae6-057c-46cf-9258-5168dcc7b0a7
# ╟─052494d1-526f-45ef-a27c-f2e7029d5b4d
# ╟─6ca184b5-133b-4c9a-a691-3a59da8663f5
# ╟─f195aae0-6e8b-41b8-a32b-4325bbf2f41b
# ╠═5fb5b59d-97a9-4292-95b1-b957d39823af
# ╟─82e023e7-e03c-4c6a-a34a-aa690d15539f
# ╟─0dacebcc-017f-4bb4-92cb-1749fd574f27
# ╟─1c133ded-ff67-487a-85a9-ec24084f0dc4
# ╟─35e5e538-3f9d-49e1-808d-8d353e02e21c
# ╟─aa71d89e-2190-4b39-bf1d-c7a680c3d1d9
# ╟─257585f4-9231-49ab-a787-ae8fa01ceb6a
# ╟─d93eeb22-e1bb-478c-a97c-7327f558800c
# ╟─e43f9ff7-7074-4129-a579-121fccd2f570
# ╠═9f631ba6-bebb-46e8-bd5c-596f5629ff27
# ╟─127b4a43-9537-42f5-a30f-10858012ceef
# ╟─dd6bec64-7b56-47a4-ade2-790008285aaf
# ╟─c02b530f-b30d-4abd-8498-8b2689078cb7
# ╟─48cde34c-dfa4-4eb6-9b2d-901a2363bacb
# ╟─ef2b1405-63c3-4f52-aaef-4adf9ae37ddb
# ╟─54d48c4f-5624-47e4-b752-fd161857fac0
# ╟─95f05467-2968-4fb5-b79d-353e11b5eaa9
# ╟─2ccbfaf8-6283-49fd-9abc-f1127eddb11f
# ╟─1f7a9428-8c7b-48c1-928c-b4f33909a79a
# ╟─aa9be589-48c6-41ee-bb29-39114d6582c5
# ╟─88426bfc-4043-43c2-92a0-44f7a22fee4e
# ╟─ff3c39a4-a447-4b34-97f1-1ff2a7736e69
# ╟─cb9f10ce-f0d2-49da-b835-8e663545c116
# ╟─e750ee9e-cce0-42fd-9e13-a01750d9fbde
# ╟─142ba9b5-9239-40ac-8d35-bc1c558d3016
# ╟─6bf2a94d-ff43-41b4-a867-c148cf8abbd6
# ╟─7acb55f0-1632-446f-b173-f15e6f56f67a
# ╟─c82157cd-b834-4871-8cfa-d5f32af7182c
# ╟─a6dacc52-83ac-4f80-96de-ac1521246592
# ╟─efd5daf8-5040-4c3e-96a9-ea71268ec43a
# ╟─f41a230e-8425-441c-a701-647fd3ede72f
# ╟─15837e7d-b35d-4e11-bb4c-e24047771e46
# ╟─eb7eeff2-5314-44d5-b10a-e07b1c68a600
# ╟─a3d63f98-6f4f-4792-b6f6-1299ca1a5a02
# ╟─f1d5822e-d077-40d1-b910-7098fd133a7d
# ╟─8f84a0e2-0faa-4e24-abd0-ff82fb910557
# ╟─17b42861-dd59-4a3c-a6a3-6ebe6806a596
# ╟─f5732129-f37e-4300-b31d-8d1ec72ea697
# ╟─5441a648-a176-4dc1-8f58-dd24e6b459c6
# ╟─67bbab74-fb4e-46a8-8ed4-12e8ef298059
# ╟─7764737a-d571-46b2-9617-d89d0d10cb04
# ╟─6f8ce07d-7128-498a-a289-170b64e446af
# ╟─edcb0089-33b9-44aa-bea6-b633840dbe66
# ╟─cf860bf9-fca1-4fd3-a1f8-465cf19f0273
# ╟─4280d95b-88a5-4a4e-822a-fe8fc4e2b0b0
# ╟─0ed9822d-3945-4525-adea-6042c2c4acfd
# ╟─8625a145-6ceb-43fb-bdb8-1368e0ad62e7
# ╟─67ea193d-ace4-4264-b20f-35724ff5db46
# ╟─30db840c-6393-4922-82c2-a147c09b24e4
# ╟─9ad4900b-c1c8-4fab-8de8-005caa983b53
# ╟─0d7c01fa-5e3f-482f-aa4b-1824ef00a275
# ╟─b6f990b7-81ec-45e6-9662-a85da2287f23
# ╟─816cb29e-b8a6-4579-a755-93709b22e967
# ╟─572cf4cd-e2ef-4469-a86d-bb757048f527
# ╟─cd4addf0-2a3a-4aa2-82ad-f17f726e143d
# ╟─24685f7d-b140-44a9-b9cf-d6d9a9471142
# ╟─faf9e517-1f09-4ea9-8e6e-e5f13863e1a9
# ╟─58e1580d-b594-4c90-aa44-897399e323f0
# ╟─e753b964-95a0-4614-830d-0b1b4ebef5c6
# ╟─22a4e90d-17c2-4187-8d8a-aac8f3601975
# ╠═4466300d-4657-4b34-bd3c-272a19550f9b
# ╟─61e4180b-f082-4174-b06b-6ed2d1759f4d
# ╠═25ee97bc-8efe-4471-9199-50291361d377
# ╠═05f46635-15bd-4132-b6cb-cee94e6f2ee3
# ╠═662be44c-48d9-4cd6-9993-984e67e74869
# ╠═2e65a895-6053-4768-9ba3-3115fef2ddac
# ╟─021eb2b3-6243-492c-8140-702a05d05350
# ╠═ae6d3d87-97ba-4d32-a560-004a94c6d506
# ╠═1dd9a5c4-70c4-4435-a679-9f78716a8c7f
# ╠═d45ebb00-2292-4c29-a13f-29f080774585
# ╠═10884f5a-1af7-4448-9d25-3efc22aff04e
# ╠═80e4ba7d-9493-41c0-a177-b06c988a86bc
# ╠═01f26215-273a-4bad-b442-578038df4724
# ╠═ee748dde-2466-4b46-816f-938abb217f09
# ╠═e22cd35d-d63f-4740-89b0-2e72b267af9e
# ╠═d8e73ce1-4444-4270-87be-5e5a619b49b1
# ╠═a8c6c34c-b6d6-4f77-b6cd-ebe9f8776271
# ╟─90c88c7d-ce7d-4db8-baeb-2cd3698b239a
# ╠═e62a5f89-47eb-4841-bde8-2b018eae0fba
# ╠═d24d621e-2350-429b-8def-dade2b04517f
# ╠═dbec0cc5-61a2-4dec-bdce-663862ceb90f
# ╠═38e6b553-63be-462a-88c4-eea9854ba73e
# ╟─b25185c9-3e1b-4cc9-9b24-54ac7492a0a1
