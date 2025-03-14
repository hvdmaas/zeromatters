extensions [ palette ]

globals [
  sum-of-spins ;; measures magnetization
  sum-of-zeros ;; counts zero spins
  mean-degree
  step ;; interval between spin values
  p-101 ;; proportion of non contrained nodes (p-101 = 1 - p-10 - p01)
  dir ;; direction of change in parameter for bifurcation plot
]

turtles-own [
  spin ;; actual spin value of a node
  spins ;; possible node values per node
]


to setup
  clear-all
  reset-ticks

  if scenario != "none" [scenarios] ;; scenerios set in interfacem linked to paper

  if not random_seed [random-seed 1] ;; if not random-seed use same initial states

  ;; Set network type including setup of turtles (nodes)
  if network = "Lattice" [Lattice]
  if network = "Erdős–Rényi" [Erdős-Rényi]
  if network = "Prefential-attachment" [Prefential-attachment]


 ;; Compute step size (interval between node value) and set spins
   set step  2 / (nr-of-spin-values - 1)
   ask turtles [set spins (range -1 (1 + step) step)]

  ;; special settings for Ising
   if model = "Ising" and nr-of-spin-values mod 2 != 0
    [ask turtles [set spins remove 0 spins]
     set nr-of-spin-values nr-of-spin-values - 1
     set alpha 0
    ]

  setup-plots

 ;; include constrained (unipolar) nodes (for depression model)
 ;; p01 and p-10 are proportions of 01 and -10 nodes
  let pp p01 + p-10
  if pp > 0 and model = "Ising" [error "Ising model not allowed when p01 or p-10 > 0"]
  if pp > 1
  [
    set p01 p01 / pp
    set p-10 p-10 / pp
  ]
  set pp p01 + p-10
  set p-101 1 - pp

  ask turtles ;; set spins and spin if p01 + p-10 > 0
  [
    let rr random-float 1
    ifelse rr < p01
    [set spins filter [x -> x >= 0] spins]
    [if rr < p01 + p-10
      [set spins filter [x -> x <= 0] spins]]

    set spin one-of spins ;; inititialize spin
    reshape ;; reshape unipolar nodes with smylies)
    recolor ;; color spins based on spin
   ]

;; for PANAS scenario create two stochastic blocks
if scenario  = "N: PANAS case: 01 and -101 nodes mixed"
  [
  ask turtles [
  ask my-links [
    let other-turtle other-end
        if ([spins = [-1 0 1]] of other-turtle and [spins = [0 1]] of myself) or
           ([spins = [0 1]] of other-turtle and [spins = [-1 0 1]] of myself)
        [if random-float 1 < 1 [die]]; deletes the link ]
    ]
  ]

    repeat 100 [  layout-spring turtles links 0.1 3 3 ]
  ]


  ;; set rho of BEG model to zero for no BEG models
  if model != "BEG" [set rho 0]

  set sum-of-spins sum [ spin ] of turtles
  set sum-of-zeros sum [ convert-boolean (spin = 0)] of turtles
  set mean-degree mean [count link-neighbors] of turtles

  set dir 1
  set scenario  "none" ;; reset scenario so that it can be changed

end

to layout
  layout-spring turtles links 0.1 30 3 ;; put linked nodes together
end


to go
  ;; update 1000 patches at a time
  repeat 1000 [
    ask one-of turtles [ update ]
  ]
  tick-advance 1000  ;; use `tick-advance`, as we are updating 1000 patches at a time
  update-plots       ;; unlike `tick`, `tick-advance` doesn't update the plots, so we need to do so explicitly
  if time-varying != "none" [time-vary]
end

;; for bifurcation plots
to time-vary
  if time-varying = "tau"
  [
    if tau > max-tv [set dir  -1]
    if tau < min-tv [set dir  1]
    set tau tau + dir * 10 ^ (-1 * c)
  ]

 if time-varying = "alpha"
  [
    if alpha > max-tv [set dir  -1]
    if alpha < min-tv [set dir  1]
    set alpha alpha + dir * 10 ^ (-1 * c)
   ]

 if time-varying = "temperature"
  [
    if temperature > max-tv [set dir  -1]
    if temperature < min-tv [set dir  1]
    set temperature temperature + dir * 10 ^ (-1 * c)
   ]


 if time-varying = "rho"
  [
    if rho > max-tv [set dir  -1]
    if rho < min-tv [set dir  1]
    set rho rho + dir * 10 ^ (-1 * c)
   ]

 if time-varying = "sigma"
  [
    if sigma > max-tv [set dir  -1]
    if sigma < min-tv [set dir  1]
    set sigma sigma + dir * 10 ^ (-1 * c)
   ]
end


to update

 ;; MCMC or Glauber updates
  let new-spin one-of remove spin spins

  let m sigma * sum [ spin ] of link-neighbors
  let q sum [ spin ^ 2] of link-neighbors
  let n count link-neighbors
  let Espin spin - new-spin
  let Espin2 spin ^ 2 - new-spin ^ 2

  let Ediff 0
  ifelse model = "Ubics"
  [  set Ediff Espin * m + tau * Espin + alpha ^ 2 * Espin2 * (q - n) ]
  [  set Ediff Espin * m + tau * Espin - alpha ^ 2 * Espin2 + rho * Espin2 * q ]

  if (Ediff <= 0) or (temperature > 0 and (random-float 1.0 < exp ((- Ediff) / temperature))) [
 ; let p min list 1 exp (- Ediff / temperature)
 ; if  random-float 1.0 < p [

;; update statistics
    let old-spin spin
    set spin new-spin
    set sum-of-spins (sum-of-spins + spin - old-spin)
    set sum-of-zeros sum-of-zeros - convert-boolean (old-spin = 0) + convert-boolean (spin = 0)]
    recolor
end

to reshape  ;; for 01 and -10 nodes (scenario PANAS)
  if max spins <= 0 [ set shape "face happy" ]
  if min spins >= 0 [ set shape  "face sad" ]
end

;; color the patches according to their spin
to recolor  ;; patch procedure
   let p-spin ((nr-of-spin-values - 1 ) / 2  * spin) + (nr-of-spin-values - 1) / 2
  set color palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 9 p-spin  0 (nr-of-spin-values - 1)

end

;; a measure of magnetization, the average of the spins
to-report magnetization
  ifelse count turtles > 0
  [report sum-of-spins / count turtles]
  [report 0]
end

;to-report variance-of-spins
;  report (sum-of-spins-squares - (sum-of-spins ^ 2) / nr-nodes) / (nr-nodes - 1)
;end

to-report zeros
  ifelse count turtles > 0
  [report sum-of-zeros / count turtles]
  [report 0]
end

to-report convert-boolean [boolean-value]
  ifelse (boolean-value) [ report 1 ] [ report 0 ]
end


to-report flatten [nested-list]
  report reduce sentence nested-list
end


to Erdős-Rényi
  resize-world -50 50 -50 50
  if nr-nodes > 500 [set nr-nodes 500]
   create-turtles nr-nodes [ set color blue
    set size 2
  set shape "circle"]
  ;repeat 10 [ask turtles [ create-link-with one-of other turtles ]]
  ;; lay it out so links are not overlapping

    ask turtles [
    ;; we use "self > myself" here so that each pair of turtles
    ;; is only considered once
    create-links-with turtles with [self > myself and random-float 1.0 < p-connect]
  ]
  repeat 100 [ layout]

  ;; leave space around the edges
  ask turtles [ setxy 1 * xcor 1 * ycor ]
  ;; put some "walker" turtles on the network

end

to Lattice
  resize-world -50 50 -50 50
  set nr-nodes 99 * 99
  let grid-size sqrt nr-nodes
  let spacing 1
  let start-x -1 * (grid-size - 1) * spacing / 2
  let start-y -1 * (grid-size - 1) * spacing / 2

  ; Create turtles in a grid
  foreach range grid-size [x ->
    foreach range grid-size [y ->
      create-turtles 1 [
        setxy (start-x + x * spacing) (start-y + y * spacing)
        set color blue
        set shape "square"
        set size 2

      ]
    ]
  ]

  ; Create links for lattice structure
  ask turtles [
    ; Link to the right neighbor
    let right-neighbor one-of turtles-on patch-at spacing 0
    if right-neighbor != nobody [ create-link-with right-neighbor ]

    ; Link to the top neighbor
    let top-neighbor one-of turtles-on patch-at 0 spacing
    if top-neighbor != nobody [ create-link-with top-neighbor ]
  ]
end

to Prefential-attachment
  if nr-nodes > 500 [set nr-nodes 500]

  set-default-shape turtles "circle"

  ;; create a random network
  create-turtles nr-nodes [
    set color blue
    set size 2
    ]
  repeat 1 [ask turtles [ create-link-with one-of other turtles ]]
  ;; lay it out so links are not overlapping
  repeat 500 [ layout-spring turtles links 0.4 4 1 ]
  ;; leave space around the edges
  ask turtles [ setxy 0.95 * xcor 0.95 * ycor ]
  ;; put some "walker" turtles on the network
end

to scenarios

 if scenario = "A: Self-organizing spatial patterns in Ising lattice model"
  [
  set network "Lattice"
  set model "Ising"
  set nr-nodes 200
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 2
  set p01 0
  set p-10 0
  set temperature .1
  set tau 0
  set alpha 0
  set rho 0
  set time-varying "none"
  set max-TV 0
  set min-TV 0
  set c 1
  ]


  if scenario = "B: Hysteresis (tau), Ising model"
  [
  set network "Erdős–Rényi"
  set model "Ising"
  set nr-nodes 50
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 2
  set p01 0
  set p-10 0
  set temperature 1
  set tau -3
  set alpha 0
  set rho 0
  set time-varying "tau"
  set max-TV 3
  set min-TV -3
  set c 2.5
  ]

 if scenario = "C: Pichtfork (temperature), Ising model"
  [
  set network "Erdős–Rényi"
  set model "Ising"
  set nr-nodes 100
  set p-connect .05
  set sigma 1
  set nr-of-spin-values 2
  set p01 0
  set p-10 0
  set temperature 0
  set tau 0
  set alpha 0
  set rho 0
  set time-varying "temperature"
  set max-TV 10
  set min-TV 0
  set c 2
  ]


  if scenario = "D: Cramer model of depression (01 Ising model)"
 [
  random-seed 1
  set network "Erdős–Rényi"
  set model "BC" ; but -1 values removed
  set nr-nodes 25
  set p-connect .2
  set sigma 1
  set nr-of-spin-values 3
  set p01 1
  set p-10 0
  set temperature .25
  set tau -3
  set alpha 0
  set rho 0
  set time-varying "tau"
  set max-TV 0
  set min-TV -5
  set c 2.5
  ]

 if scenario = "E: Hysteresis in Ising model with 100 spin values per node"
  [
  set network "Erdős–Rényi"
  set model "Ising"
  set nr-nodes 50
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 101
  set p01 0
  set p-10 0
  set temperature 1
  set tau -1.5
  set alpha 0
  set rho 0
  set time-varying "tau"
  set max-TV 1.5
  set min-TV -1.5
  set c 3
  ]

 if scenario = "F: Tricritical transition (figure 2)"
  [
  random-seed 1
  set network "Erdős–Rényi"
  set model "BC"
  set nr-nodes 100
  set p-connect .05
  set sigma 1
  set nr-of-spin-values 3
  set p01 0
  set p-10 0
  set temperature .5
  set tau 0
  set alpha 1
  set rho 0
  set time-varying "alpha"
  set max-TV 2.4
  set min-TV 1
  set c 2.5
  ]

 if scenario = "G: Stochastic resonance in BC at low alpha"
  [
  random-seed 1
  set network "Erdős–Rényi"
  set model "BC"
  set nr-nodes 100
  set p-connect .05
  set sigma 1
  set nr-of-spin-values 3
  set p01 0
  set p-10 0
  set temperature 2.5
  set tau 0
  set alpha 0
  set rho 0
  set time-varying "none"
  set max-TV 0
  set min-TV 0
  set c 1
  ]


  if scenario = "H: Tricritical behavior in BC"
  [
  random-seed 1
  set network "Erdős–Rényi"
  set model "BC"
  set nr-nodes 100
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 3
  set p01 0
  set p-10 0
  set temperature 2
  set tau 0
  set alpha 2.27
  set rho 0
  set time-varying "none"
  set max-TV 0
  set min-TV 0
  set c 1
  ]

  if scenario = "I: Fig 4: Double or pinched hysteresis"
  [
  set network "Erdős–Rényi"
  set model "BC"
  set nr-nodes 50
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 3
  set p01 0
  set p-10 0
  set temperature .35
  set tau -2
  set alpha 1.7
  set rho 0
  set time-varying "tau"
  set max-TV 2
  set min-TV -2
  set c 2.5
  ]

  if scenario = "J: Tricritical behavior in Ubics in Preferential attachment model"
  [
  random-seed 1
  set network "Prefential-attachment"
  set model "Ubics"
  set nr-nodes 200
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 3
  set p01 0
  set p-10 0
  set temperature .350
  set tau 0
  set alpha 1
  set rho 0
  set time-varying "none"
  set max-TV 0
  set min-TV 0
  set c 1
  ]

  if scenario = "K: Trimodal distribution of m in Ubics lattice model"
  [
  set network "Lattice"
  set model "Ubics"
  set nr-nodes 200
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 101
  set p01 0
  set p-10 0
  set temperature .1
  set tau 0
  set alpha .82
  set rho 0
  set time-varying "none"
  set max-TV 0
  set min-TV 0
  set c 1
  ]

  if scenario = "L: Mixed islands in the BEG lattice model"
   [
  set network "Lattice"
  set model "BEG"
  set nr-nodes 200
  set p-connect .1
  set sigma 0
  set nr-of-spin-values 51
  set p01 0
  set p-10 0
  set temperature .1
  set tau 0
  set alpha 1.25
  set rho 1
  set time-varying "none"
  set max-TV 0
  set min-TV 0
  set c 1
  ]



  if scenario = "M: Rho bifurcation"
  [
  set network "Erdős–Rényi"
  set model "BEG"
  set nr-nodes 50
  set p-connect .1
  set sigma 1
  set nr-of-spin-values 3
  set p01 0
  set p-10 0
  set temperature 2
  set tau 0
  set alpha 2
  set rho 0
  set time-varying "rho"
  set max-TV 2
  set min-TV 0
  set c 3
  ]



 if scenario = "N: PANAS case: 01 and -101 nodes mixed"
 [
  random-seed 1
  set network "Erdős–Rényi"
  set model "Ubics"
  set nr-nodes 100
  set p-connect .06
  set sigma 1
  set nr-of-spin-values 3
  set p01 0.5
  set p-10 0
  set temperature .75
  set tau -1
  set alpha 1
  set rho 0
  set time-varying "tau"
  set max-TV 1
  set min-TV -1
  set c 2.5
  ]



  if scenario = "Asymmetric hysteresis in BC"
  [
  set network "Erdős–Rényi"
  set model "BC"
  set nr-nodes 20
  set p-connect 1
  set sigma 1
  set nr-of-spin-values 7
  set p01 0
  set p-10 0
  set temperature .1
  set tau 0
  set alpha 4
  set rho 0
  set time-varying "alpha"
  set max-TV 5
  set min-TV 1
  set c 3
  ]




end
@#$#@#$#@
GRAPHICS-WINDOW
315
10
893
589
-1
-1
5.644
1
10
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
attempted flips
1000.0

BUTTON
160
10
305
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
11
470
221
503
temperature
temperature
0
5
1.0
0.01
1
NIL
HORIZONTAL

PLOT
910
15
1445
265
Magnetization
time
average spin
0.0
20.0
-1.1
1.1
true
true
"" ""
PENS
"m (magn.)" 1.0 0 -13345367 true "" "plotxy ticks magnetization"
"n ('zeros')" 1.0 0 -8732573 true "" "plotxy ticks zeros \n;plotxy ticks variance-of-spins "

BUTTON
10
10
155
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
11
505
306
538
tau
tau
-12
12
-0.5050000000000546
.01
1
NIL
HORIZONTAL

SLIDER
15
150
310
183
nr-nodes
nr-nodes
0
500
50.0
1
1
NIL
HORIZONTAL

SLIDER
11
540
306
573
alpha
alpha
0
8
8.0
.01
1
NIL
HORIZONTAL

PLOT
1190
280
1430
450
Histogram
magn.
NIL
-1.0
1.0
0.0
10.0
true
false
"set-plot-x-range -1 - .05 1 + .05\n;set-plot-x-range min-spin - .5 max-spin + .5\n;set-plot-y-range 0 count turtles\nprint nr-of-spin-values\nset-histogram-num-bars nr-of-spin-values" ""
PENS
"default" 1.0 1 -5825686 true "" "histogram [spin] of turtles\n;histogram [round max-spin * spin] of turtles"

SLIDER
15
210
190
243
p-connect
p-connect
0
1
0.1
.001
1
NIL
HORIZONTAL

PLOT
940
280
1180
450
m,n
NIL
NIL
-1.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 2 -16777216 true "" "plotxy magnetization zeros"

CHOOSER
15
100
197
145
network
network
"Lattice" "Erdős–Rényi" "Prefential-attachment"
1

SWITCH
15
60
157
93
Random_seed
Random_seed
0
1
-1000

CHOOSER
175
100
313
145
model
model
"Ubics" "BC" "Ising" "BEG"
2

SLIDER
11
580
306
613
rho
rho
0
1
0.0
.01
1
NIL
HORIZONTAL

SLIDER
14
310
289
343
nr-of-spin-values
nr-of-spin-values
3
1001
100.0
2
1
NIL
HORIZONTAL

SLIDER
15
404
226
437
p01
p01
0
1
0.0
.01
1
NIL
HORIZONTAL

MONITOR
247
368
304
413
  p-101
p-101
3
1
11

MONITOR
215
230
307
275
mean degree
mean-degree
0
1
11

TEXTBOX
30
352
239
379
prop. of [-1..0] and [0..1] spins
11
0.0
1

MONITOR
236
455
293
500
Beta
1 / temperature
2
1
11

SLIDER
15
265
187
298
sigma
sigma
0
1
1.0
.01
1
NIL
HORIZONTAL

TEXTBOX
15
250
165
268
connection strength:
11
0.0
1

TEXTBOX
20
195
170
213
if \"Erdős–Rényi\":
11
0.0
1

TEXTBOX
31
440
208
467
Model parameters. Can be varied during simulation:
11
0.0
1

TEXTBOX
325
595
475
613
red to blue = [-1,1]
11
0.0
1

CHOOSER
30
650
168
695
time-varying
time-varying
"none" "tau" "alpha" "temperature" "rho" "sigma"
1

SLIDER
25
800
209
833
c
c
1
10
3.0
.1
1
NIL
HORIZONTAL

PLOT
319
623
898
828
bifurcation
time-varying
m
0.0
0.0
-1.0
1.0
true
false
"" ""
PENS
"default" 1.0 2 -2674135 true "" "if time-varying  = \"tau\" and dir = 1 [plotxy tau magnetization]\nif time-varying  = \"alpha\" and dir = 1 [plotxy alpha magnetization]\nif time-varying  = \"temperature\" and dir = 1 [plotxy temperature magnetization]\nif time-varying  = \"rho\" and dir = 1 [plotxy rho magnetization]\nif time-varying  = \"sigma\" and dir = 1 [plotxy sigma magnetization]"
"pen-1" 1.0 2 -10649926 true "" "if time-varying  = \"tau\" and dir = -1 [plotxy tau magnetization]\nif time-varying  = \"alpha\" and dir = -1 [plotxy alpha magnetization]\nif time-varying  = \"temperature\" and dir = -1 [plotxy temperature magnetization]\nif time-varying  = \"rho\" and dir = -1 [plotxy rho magnetization]\nif time-varying  = \"sigma\" and dir = -1 [plotxy sigma magnetization]"

INPUTBOX
113
712
185
772
max-TV
1.5
1
0
Number

INPUTBOX
25
712
99
772
min-TV
-1.5
1
0
Number

TEXTBOX
30
628
218
657
bifurcation plots
11
0.0
1

TEXTBOX
43
780
231
803
change = 10^-c
11
0.0
1

CHOOSER
939
612
1373
657
scenario
scenario
"none" "A: Self-organizing spatial patterns in Ising lattice model" "B: Hysteresis (tau), Ising model" "C: Pichtfork (temperature), Ising model" "D: Cramer model of depression (01 Ising model)" "E: Hysteresis in Ising model with 100 spin values per node" "F: Tricritical transition (figure 2)" "G: Stochastic resonance in BC at low alpha" "H: Tricritical behavior in BC" "I: Fig 4: Double or pinched hysteresis" "J: Tricritical behavior in Ubics in Preferential attachment model" "K: Trimodal distribution of m in Ubics lattice model" "L: Mixed islands in the BEG lattice model" "M: Rho bifurcation" "N: PANAS case: 01 and -101 nodes mixed" "Asymmetric hysteresis in BC"
0

TEXTBOX
943
556
1415
602
Here you can select pre-defined scenarios that illustrate properties of the models, discussed in the paper:
13
14.0
1

SLIDER
15
368
227
401
p-10
p-10
0
1
0.0
.01
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

This is a simulation tool for the paper "The Statistical Physics of Psychological Networks: Zero matters" of van der Maas et al. (2025).
It allows readers new to this type of modeling to visually verify the properties of several spin models. The program incorporates the Ising, Blume Capel (BC), Ubics, Blume-Emery-Griffiths (BEG) model, different network types (Lattice, Erdős–Rényi, prefential-attachment), and allows one to vary network size, connectivity and connection strength. In addition, spin values can take on any number between -1 and 1 (the general-spin case). Also, mixtures of node types (some being only negative or positive) can be investigated. 


## HOW IT WORKS

We represent the possible spin states with numbers between +1 or -1. Spins of +1 are shown in blue, spins of -1 in red. Intermedia colors represent intermediate values. Light-yellow is the zero state.

A spin decides whether to "flip" to its opposite as follows. The spins are seeking a low energy state, so a spin tend to flip if flipping would decrease its energy. But the spins sometimes also flip into a higher energy state if the temperatue is high. We calculate the exact probability of flipping using the Metropolis algorithm, which works as follows. Call the potential gain in energy Ediff. Then the probability of flipping is:

e<sup>-Ediff / temperature</sup>.

The gist of this formula is that as the temperature increases, flipping to a higher energy state becomes increasingly likely, but as the energy to be gained by flipping increases, the likelihood of flipping decreases. You could use a different formula with the same gist, but the Metropolis algorithm is most commonly used.

To run the model, we repeatedly pick a single random spin and give it the chance to flip.

## HOW TO USE IT

The main model parameters a,T (1⁄β),τ, and ρ  can be varied while nodes are continuously updated. This way the user can move over the phase diagram (figure 3 in the paper) and construct bifurcation plots. Specific scenarios are readily accessible and will be referenced throughout the paper.

Then press GO to watch the model run.

The magnetization of the system is the average (mean) of all the spins. The MAGNETIZATION monitor and plot show you the current magnetization and how it has varied so far over time.

## THINGS TO NOTICE

See scenarios referenced in the paper

## THINGS TO TRY

Start with any scenario and adapt using the sliders


## HOW TO CITE

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment-40nodes_t=2" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>spinofturtles</metric>
    <steppedValueSet variable="alpha" first="0.8" step="0.01" last="1.1"/>
    <steppedValueSet variable="tau" first="-0.35" step="0.02" last="0.35"/>
    <enumeratedValueSet variable="initial-zeros">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-connect">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nr-nodes">
      <value value="40"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="1run" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="20000"/>
    <metric>magnetization</metric>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-zeros">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-connect">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nr-nodes">
      <value value="40"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_12nodesb" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="50"/>
    <metric>spinofturtles</metric>
    <steppedValueSet variable="alpha" first="0.6" step="0.02" last="1.4"/>
    <steppedValueSet variable="tau" first="-0.5" step="0.025" last="0.5"/>
    <enumeratedValueSet variable="initial-zeros">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-connect">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nr-nodes">
      <value value="20"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
