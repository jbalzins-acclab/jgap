# AI usage in JGAP development

This document is created to add some transparency regrding the AI usage throughout the development of JGAP.
Early prototypes (should be available through Git history) we created with a rather minimalistic AI usage,
namely ChatGPT + Claude were used primarily as the replacement of what used to be Google (+stack overflow...) 
workflow, Clion's lightweight autocomplete and early version of JetBrains Junie drafted the logging.
Later on, the project went through a lot of, first just refactoring, later almost complete overhauls.

At some point in (approximately) April 2026, the auther got familiarized with agentic AI, and Google's Gemini started to be used actively,
however, with extreme caution (adjusting the code only by small pieces, manual compilation and testing).
Later on, once the core architecture seemed to be done, Calude agent helped more extensively with serialization;
some basic ideas and patterns remained from older versions (Registry), but serializers were written completely by Claude.

Such a workflow felt unreliable, so proceeding adjustments to core architecture were done 
with Gemini again asking to edit small, although at that point larger pieces.
The overall structure/architecture of the core, main logic and workflow are the author's own.
Gemini did, however, help with brainstorming and feedback, and it's largest contribution was the 
suggestion of the algorithm for solving least-squares block-incrementally (it took quite a long while
for it to propose it, and a bunch of methods like randomized SVD and iterative normal equation fitting were 
attempted first with the goal of tried to optimize memory usage; the algorithm itself is not due AI, 
but is relatively old and called: 
"banded sequential Householder triangularization" - see C. L. Lawson and R. J. Hanson. In: Solving Least Squares Problems. Society for
Industrial and Applied Mathematics, Jan. 1995, pp. 10–17, 207–232.),
while element-incremental fit was completely author's idea based on old observation about high number of zeroed blocks in the covariance matrices.

Early unit tests were hand crafted, however, around the time Serialization was reimplemented, AI took over largely, with author just veryfing new tests.
Compilation was AI guided, with author focusing only on the end-goal, steering AI into something that would work.
Similarly with Python interface and publishing details.
Biggest fraction of the DOCS top level docs is the gibberish in the author's head translated into understandable ReadMe's,
while compilation, installation and convinience stuff had more AI authorship to it.
In code comments in src/jgap are primarily author's while AI used to describe e.g. usage in all other files.
Example, validation, and key unit test logic are of author's design.