//Gauss Elimination by Asmaa Ali - ID: 1910069

using namespace std;

class GaussElimination
{
    public:
        /**
        @param float **a  the augmented matrix
        @param  n  array size
        */
        void printMatrix(float **a, int n);
        /**
        @param float **a  the augmented matrix
        @param  n  array size
        @param  scale the scale vector
        */
        void selectMax(float **a, int n, vector<float> &scale);
        /**
        @param float **a  the augmented matrix
        @param  n  array size
        */
        void backwardSubstitution(float **a, int n);
        /**
        @param float **a  the augmented matrix
        @param  n  array size
        @param  scale the scale vector
        @param  step the forward elimination step
        */
        void solve(float **a, int n, vector<float> &scale, int step);
        /**
        Checks solution validity
        */
        bool valid_solution = false;
    protected:

    private:
        /**
        @param float **a  the augmented matrix
        @param  n  array size
        @param  scale the scale vector
        @param  step the forward elimination step
        */
        void scaledPartialPivoting(float **a, int n, vector<float> &scale, int step);
        /**
        @param float **a  the augmented matrix
        @param  n  array size
        @param  step the forward elimination step
        */
        void forwardElimination(float **a, int n, int step);
};
