Fan schedule files for REGCAP simulations (RIVEC and COMMISSIONING) WJNT 2011

3 different files:

schedule1*.txt for Prototype B (1,200 sq ft house with 2 bathrooms and 4 occupants)
schedule2*.txt for Prototype C (2,100 sq ft house with 3 bathrooms and 4 occupants)
schedule3*.txt for Prototype D (2,400 sq ft house with 3 bathrooms and 5 occupants)

* denotes a,b or c. These are the same files with different names. Used to overcome an issue with the Qbasic version of REGCAP where running different instances of REGCAP using the same schedule input file caused issues.

schedule?.csv (where ? = 1,2,3) is the comma separated version of the same file. The C++ version of REGCAP reads in the tab delimited (.txt) versions.

The files are a list of 0's and 1's for each minutes of the year. 0 = fan OFF. 1 = fan ON.

Column 1 = dryer fan
Column 2 = kitchen fan (range hood)
Column 3 = bathroom fan 1
Column 4 = bathroom fan 2
Column 5 = bathroom fan 3

Update 2/10/16 LIR
REGCAP now reads in one file since c++ does not have lock files for read.
File names are now sched1.txt, sched2.txt, and sched3.txt.
Added 5th column to sched1.txt with all 0's so that same routine can read it in.
