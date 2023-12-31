# Perseverative Thought
Code for dynamic and control theory analysis of brain activity during worry and rumination in depression and anxiety patients

Adapted from code in:

Cornblath EJ, Ashourvan A, Kim JZ, Betzel RF, Ciric R, Adebimpe A, Baum GL, He X, Ruparel K, Moore TM, Gur RC. Temporal sequences of brain activity at rest are constrained by white matter structure and modulated by cognitive demands. Communications biology. 2020 May 22;3(1):261.

https://github.com/ejcorn/brain_states

Singleton SP, Luppi AI, Carhart-Harris RL, Cruzat J, Roseman L, Nutt DJ, Deco G, Kringelbach ML, Stamatakis EA, Kuceyeski A. Receptor-informed network control theory links LSD and psilocybin to a flattening of the brain’s control energy landscape. Nature communications. 2022 Oct 3;13(1):5812.

https://github.com/singlesp/energy_landscape

## Main analysis

Code > Analysis > 

1. Preprocess and run k-means clustering data A0_preprocess_brainStates.m
2. Plot clustering result A0_plot_brainStates.m
	- Finding partition for rest of analysis:
		- repeatkmeans_sps.m
		- elbow_sps.m
		- ami_calc.m
3. Test reliability of clustering across individuals and time A0_reliability_brainStates.m
4. Group and within-person analysis with worry and rumination combined Code > Analysis > worryANDruminate > 
	- plot_dwellTime.m
	- plot_fractionalOccupancy.m
	- plot_transitionProbability.m
5. Group and within-person analysis with worry and rumination separated Code > Analysis > worryXORruminate >
	- A0_preprocess_brainStates.m
	- plot_dwellTime.m
	- plot_fractionalOccupancy.m
	- plot_transitionProbability.m
6. Control energy analysis plot_controlEnergy.m