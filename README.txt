---------------------------------------
IS DRUMMING FRACTAL?
OLIVER MILES GORDON & DOMINIC ARRAN COY
---------------------------------------

- 2 MAIN FILES: Fin_Main.m, Fin_Analysis.m
- 1 FUNCTION: filterlims.m
- 1 DATABASE FILE: \Docs\People.xlsx
- 1 ORIGINAL TOM SAWYER FILE: \Samples\Tom Sawyer For Musician.wav

---------------
TO ADD NEW DATA
---------------
- Set up microphone in mono, at 44.1kHz
- Record musician playing in double handed, single handed, metronomic technique
	while listening to Tom Sawyer
- Recording of Tom Sawyer at \Samples\Tom Sawyer For Musician.wav
- They hear 2 bars before piece starts. On first bar, do nothing.
	On second bar, play 4 crotchets.
	Then play semiquavers
- Clip off start of recording before crotchets start
- Assign person a number, and fill in database \Docs\People.xlsx
- Save recording as .wav, with name Metronomic 1.wav, Double Handed 1.wav, Single Handed 1.wav
	(but replace number with number of musician)

--------------
Fin_Main.m
--------------
- Settings for changing are at start
- sample_start and sample_end are numbers of samples to analyse
- technique should be either "Single Handed","Double Handed","Metronomic" (without quotations)
	to analyse that particular drumming technique
- All other settings shouldn't need to be changed
- Various options to test DFA algorithm around line 151:
	- Uncommenting line 48,152/153 selects data from Rasenen paper (should reproduce results)
	- Uncommenting lines 154-156 generates noise with selectable Hurst exponent (end of
		line 155). Need to comment lines 163-166 (Should output result input)
	- Uncommenting line 157 generates data as if from perfect metronome
- Can add noise by uncommenting lines 169-170
- If wanting to do DFA on power, not time, rename delta_p as delta_t before line 178

--------------
Fin_Analysis.m
--------------
- Lines 10 and 18 need to be altered depending on analysis to be done:
	- If wanting just double handed, make these lines equal 1 and 1
	- If wanting just single handed, make these lines equal 2 and 2 (and then manually x2 output numbers)
	- If wanting just metronomic, make these lines equal 3 and 3 (and then manually x3 output numbers)
	- If wanting combined single+double handed, make these lines equal 2 and 1:2
	- If wanting combined single+double handed+metronomic, make these lines equal 3 and 1:3
	- (sorry)

--------------
filterLims.m
--------------
- Nothing needs changing here

-----------------------------------
- Thank you for your help and time
	 - Dominic & Oliver
-----------------------------------