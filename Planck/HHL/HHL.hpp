namespace Planck {

    class HHL{
        public:
            HHL(matrix A, vector b, int m);
            vector evaluate();

        private:
            matrix A; 
            vector b;
            int m;
            Circuit circuit;

    };

}