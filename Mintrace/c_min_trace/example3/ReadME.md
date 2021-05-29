I have finished the Householder and QR function in TraceSolver, which are devised from example2;
What's more, there is a bug:
If just create a default CGSolver like this:
CGSolver AAA;
then make this .cpp will be error;
but add command like:
SparseMatrix<double> Temp;
After above command, this .cpp file will work well!
