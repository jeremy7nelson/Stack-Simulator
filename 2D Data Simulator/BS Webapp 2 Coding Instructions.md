date: 2026-08-05

- Forgot to account for scattering of different models across the surface. 

# Introduction / Instructions
The following are instructions for generating a script in Javascript. The final product should be a .js type file. All parameters are to be hard-coded, such that a different person can write an interface for it, i.e. turn this into an html document. Let there be a separate section for those hard-coded parameters that are to be user-defined, and let there be comments following those assignments that indicate what form of selection should take place. i.e. an empty box to be filled, a choice from a list, or a box to be checked or unchecked. Furthermore, there should be some indication in the comments of when user-defined parameters are conditional on model selection, and when they are not, as guidance for a another person looking to design an interface.

I want to combine elements of existing work into an new script, as follows:
### Existing work: 
spr_simulator_draft8.html
- This html file is an SPR data simulator that produces single, kinetic time data series for a chosen model. 
- This program incorporates all models under discussion: the Langmuir, the Heterogeneous Analyte, the Conformational Change, the Bivalent Analyte, as well as the option to enable mass transport limitation for any of the four. We want this in the new program. 
- This program is wrapped in an HTML visualization. We do *not* want that. 

generate_data_array_tesselate2.js: 
- This script is the current iteration of the two-dimensional data simulator. At present it generates $s[k]$ in the appropriate manner, including the assigning of capacity bins and the creation of the $\text{binField}$ vector. We want these elements in the new script.
- This script incorporates a scheme for pseudo-random generation of spatial elements and noise elements where each of the two starts from a unique seed. This scheme is to be preserved. The default seeds may be used again in the new script. 
- This script incorporates an Ornstein-Uhlenbeck process for common drift that serves as the starting point for modifications described below. The same OU parameters found within generate_data_array_tesselate2.js may be used here. 
- This script assembles the final stack in a two-step method, building the noise separately from the date before combining them. That is being retired in favor of the scheme illustrated below under the "Map" heading. 
- This script only models the Langmuir dynamic. 
- This program integrates all kinetic data at unit capacity. We are not adhering to this convention. See details below, under the heading "Rmax."

"BI - 2Dimensional Data Simulator Documentation.md"
- Documentation for an earlier version of generate_data_array_tesselate2.js.
- Valuable primarily for defining which variable are defined by the user, and which variables are not. Use in conjunction with the definitions under the "definitions" heading to see find the sum total of what should be user-defined. 
### Details Regarding the Heterogeneous Ligand Model
In the case of the Heterogeneous Ligand Model:
- Adopt the conventions found in spr_simulator_draft8.html with respect to the integration of a kinetic data series both with and without the effects of Mass Transport Limitation. i.e.
	- Both $R_{\text{max},1}$ and $R_{\text{max},2}$ are to be defined by the user.
	- When it comes to default values, however, let $R_{\text{max},1}=R_{\text{max},2}=0.5 \cdot \text{RmaxD}$
### Rmax
Rather than integrate at unit capacity, we will now integrate according to some user-defined value $R_{\text{max}}$.
- We will need a new, appropriate default value for $k_{\text{tr}}$. Use $1\text{e}9$,
- We will need an appropriate default value for $R_{\text{max}}$, which in spr_simulator_draft8 is 120.
- Let the default value of $R_{\text{max}}$ be assigned to a variable, similar to ```const RmaxD = 120```, or whatever the appropriate default value is, such that:
	- In the case of jitter, instead of $\sigma = 0.10 R_{max}$, let $\sigma = 0.10 *\text{RmaxD}$ 
	- In the case of drift, let likewise $D = 3.5 * \text{RmaxD}$
- The definition if $s[k]$ does not change. The capacity bins do not change. The what must be different is that, integration is against $R_{\text{max}}\cdot \text{capacity}$, which is where capacity is either a single value, 1, or one of several values from $\text{capacityBins}$ in the case of bivalent analyte or MTL.
- # Map
Goal: to generate SPRm data according to the following scheme:
$$\text{Stack}[k][t] = s[k]\cdot \hat{Y}_{\text{bin}[k]}[t] + \text{nsaRate}[k]\cdot Y_{ns}[t] + g[k]\cdot\text{driftCommon}[t] + \sigma\cdot\text{hash}(i,j,t)$$
# Definitions:
- $m\times n$: the dimensions of the canvas, where $m$ is the number of rows, and $n$ the number of columns. By default $640 \times 480$. This is not a typo, that is 640 rows and 480 columns. User-defined. 
- $k$: index within the spatial vectors. Defined as $k = i \cdot n +j$, where $n$ is the number of pixels in the horizontal direction. 
- $\Delta t$: Parameter that defines resolution in time. Defaults to 1(s). Hard-coded.
- $\text{Stack}[k][t]$: This is the final product, a simulation of SPRm time-series data in an $m\times n$ window along discretized time vector $T$. Shape is two-dimensional: $(m \times n) \times \text{Length}(T)$. The entire stack is held in memory as a float32 array.
- $s[k]$: A scalar field, which is a vector with dimension $(1 \times m\cdot n)$, which encodes the circular regions, maximum capacities of circular regions, and variation within each circular region from a maximum at the circles edge to a minimum of $(0.15 \cdot \text{capacity})$ in the center. Every element of $s$ lies between $0.15\times \text{capacity}$ and $\text{capacity}$. Capacity for each circular region is pseudo-randomly selected during the construction of the of field.
- $\text{capacityBins}$: A short vector which stores the collection of discrete capacities that may be assigned to each circular region. Elements have values between 0 and 1. 
- $\text{binField}$: Constructed alongside $s[k]$. A vector which stores, for each pixel, an index which corresponds to an element of $\text{capacityBins}$. Used as a reference when applying an integrated kinetic data series to a pixel.   
- $Y_{ns}[t]$: Langmuir Series integrated at: 
	- $R_{ns} = 1$. 
	- $k_{d,ns} =$ 1e-4
	- $k_{a,ns} =$ 4e3 
- $Y_{\text{bin}[k]}[t]$: Kinetic data series. For models Langmuir and Conformational change there is one data series. For Heterogeneous Ligand there are two. For bivalent analyte and any model involving MTL, there will be one unique kinetic data series per capacity bin, per unique set of kinetic parameters. 
- $\hat{Y}_{b}= Y_{\text{b}}/\text{cap}_b$: Normalized kinetic data series - prevents a double application of the capacity restriction. 
- $\hat{g}[k]=|\nabla s|_k/p_{95}$: Sobel field representation of the gradient, normalized to the 95th percentile, such that the field is dimensionless and $R_{\text{max}}$ invariant. 
- $\text{nsaRate}[k]$: The weighting element that scales the contribution of NSA. Defined $\text{nsaRate}=a_n + b_n \hat{g}[k]$
	- $a_n$: User-defined parameter that governs the scale of the space-independent component of NSA. Should be greater than zero. Default at 80 - makes background around equal to jitter.
	- $b_n$: User-defined parameter that governs the scale of the spatial component of NSA. May be equal to zero to remove the spatial component of NSA. Default at 80: Should make edge effect about double the background. 
- $Y_{ns}[t]$: Langmuir Series modeling non-specific adsorption. Integrated at: 
	- $R_{ns} = 1$. (Let this remain 1, scale managed by $\text{nsaRate}$)
	- $k_{d,ns} =$ 1e-4 - try this first
	- $k_{a,ns} =$ 4e3 - try this first
- $g[k]$: Weight parameter for common drift. Will use isotropic drift for a first try. Defined $g[k]=a_d + b_d\hat{g}[k]$ 
	- $a_d$: The non spatial component of drift. Hard-coded at 1. 
	- $b_d$: User-Defined. Scales the spatial component of drift. May be zero to remove the spatial component of drift. Should lie between 0.03 and 0.3
- $\text{RmaxD}$ : Default, hard-coded value of $R_{\text{max}}$
- $\text{driftCommon[t]}$: Modeled as $D\left(1 - e^{-t/\tau}\right)+ \omega(t)$, where 
	- $D$ is the magnitude of the drift. By default $D = 3.5 * \text{RmaxD}$, where 
	- $\tau$: Exponential decay parameter, governs how quickly drift accumulates. Hard-coded at 500.
	- $\omega(t)$: Ornstein-Uhlenbeck component of drift that attenuates according to the same decay parameter. 
- $\sigma$: Scaling component of per-component jitter. By default $= 0.10 *\text{RmaxD}$
- $\text{hash}(i,j,t)$: Pseudo-random noise with a Gaussian distribution according to the Box-Muller Transform of uniformly distributed values generated by an integer hash of $(i,j,t,\text{seed})$.
# Integration and Assembly steps:
- Phase 0: Integrate. All the ODE work happens here, producing a handful of time-course vectors: $Y_{\text{bin}[k]}[t]$ for each capacity bin, $Y_{\text{ns}}[t]$ for NSA, and $\text{driftCommon}[t]$. How many $Y_{\text{bin}}$ vectors you need depends on whether MTL or bivalent is active (four, one per bin) or not (one, shared)
- Phase 1: Build per-pixel weights. This includes the scalar capacity field that incorporates circles, and the power diagram. Includes the derivation of every scalar needed by the inner loop. No time-axis component of this step. 
- Phase 2: Assemble each element of $\text{stack}[k][t]$ in a frame-major fashion. 
# Appendix
### Possible future additions or concerns
Possible Problems:
- Internal Boundaries
Possible additions:
- A provision for "circle owner"
- A different characterization of drift, as described below.
### Additional notes on Drift
##### Note Re Drift:
The drift is far more pronounced early in the injection series. The rate of drift decreases substantially as the experiment proceeds, and has largely settled down by the third or fourth injection. We have chosen to model that rate with an exponential decay function.
##### Two ways to model drift:
Isotropic: $g[k] = a + b|\nabla s|_k$  or Lateral Shift: $\text{drift}_k(t) = a\cdot \text{driftCommon}(t) + \nabla s_k \cdot \delta(t)$
In the case of lateral shift, delta is a two-component OU process representing sub-pixel stage wander. Represented as follows:

$$+ \text{gUniform}[k]·\text{driftCommon}[t] + g_i[k]·DX[t] + g_j[k]·DY[t]$$






