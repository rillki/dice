module eonium.la.mat;

/++ MATRIX module
    - struct Matrix
        - Matrix: from 2d array
        - Matrix: from size (rows, cols)
        - dup
        - resize
    - matDiagonal
+/

import std.algorithm: copy;
import std.traits: isFloatingPoint;

alias Matrixf = Matrix!float;
alias Matrixd = Matrix!double;
alias matDiagonalf = matDiagonal!float;
alias matDiagonald = matDiagonal!double;

struct Matrix(T = float) if(isFloatingPoint!T) {
    private {
        T[] data;
        size_t r;
        size_t c;
    }
    
    T[][] m;
    alias m this;
    
    /// Construct matrix from 2d array
    this(const T[][] arr) in (arr.length && arr[0].length) {
        // save matrix size
        r = arr.length;
        c = arr[0].length;

        // copy matrix data
        data.length = r * c;
        foreach(i, row; arr) {
            foreach(j, val; row) {
                data[i * r + j] = val;
            }
        }

        // set slices
        updateSlices();
    }

    /// Construct a matrix from size (rows, cols) and fill with value
    this(const size_t r, const size_t c, const T value = 0) in (r && c) {
        // save matrix size
        this.r = r;
        this.c = c;

        // create an empty matrix and fill it with value
        data.length = this.r * this.c;
        data[] = value;
        
        // set slices
        updateSlices();
    }
    
    /// Returns number of rows
    const(size_t) rows() const {
        return r;
    }

    /// Returns number of cols
    const(size_t) cols() const {
        return c;
    }

    /// Returns matrix data as 1d array
    const(T[]) array() const {
        return data;
    }

    /// Checks if matrix entries are 0
    bool isZero() const {
        bool zero = true;
        foreach(d; data) {
            if(d != 0) {
                zero = false;
                break;
            }
        }

        return zero;
    }
    
    /// Returns a copy of a matrix
    Matrix!T dup() const {
        return Matrix!T(this);
    }
    
    /// Resizes the matrix to (rows, cols)
    void resize(const size_t r, const size_t c) in (r && c) {
        // check if we need to allocate more memory
        if(r * c > this.r * this.c) {
            data.length = r * c;
        }

        // update matrix size
        this.r = r; 
        this.c = c;
        
        updateSlices();
    }

    /// Element wise increment/decrement
    ref Matrix!T opUnary(string op)() if(op == "++" || op == "--") {
        mixin(op ~ "data[];");
        return this;
    }

    /// Matrix operations (+, -, *)
    Matrix!T opBinary(string op)(const Matrix!T rhs) 
        if(op == "+" || op == "-" || op == "*") 
        in(((op == "+" || op == "-") && rows == rhs.rows && cols == rhs.cols) || (op == "*" && rows == rhs.cols)) 
    {
        auto result = this.dup();

        if(op == "+" || op == "-") {
            foreach(i, ref row; result) {
                foreach(j, ref entry; row) {
                    mixin("entry " ~ op ~ "= rhs[i][j];");
                }
            }    
        } else {
            //...
        }

        return result;
    }
    
    /// Element wise value operations (+, -, *, /)
    ref Matrix!T opOpAssign(string op)(const T value) if (op == "+" || op == "-" || op == "*" || op == "/") {
        mixin("data[] " ~ op ~ "= value;");
        return this;
    }

    /// Matrix operations (+=, -=)
    ref Matrix!T opOpAssign(string op)(const Matrix!T mat) if (op == "+" || op == "-") {
        mixin("data[] " ~ op ~ "= mat.array[];");
        return this;
    }

    private {
        /// Updates slice pointers
        void updateSlices() {
            m.length = this.r;
            foreach(i, ref row; m) {
                immutable s = i * this.c;
                row = data[s..(s + this.c)];
            }   
        }
    }
}

/++
Creates a diagonal matrix with custom value filled along the diagonal (default: 1)

Params:
    rows = number of rows
    cols = number of cols
    value = diagonal fill value

Returns: Matrix!T
+/
Matrix!T matDiagonal(T = float)(const size_t r, const size_t c, const T value = 1) if(isFloatingPoint!T) in (r && c) {
    auto mat = Matrix!T(r, c, 0);
    for(size_t i = 0; i < (r < c ? r : c); i++) {
        mat[i][i] = value;
    }

    return mat;
}

/* --------------------- UNITTESTS --------------------- */

unittest {
    import std;
    
    // from 2d array
    float[][] f = [[-15, 2], [3, 4]];
    auto mat = Matrix!float(f);
    auto matf = Matrixf(f);
    auto matd = Matrixd(f.to!(double[][]));
    assert(matf == f);

    // from size
    auto mat0 = Matrixf(2, 2);
    assert(mat0 == [[0, 0], [0, 0]].to!(float[][]));

    auto mat1 = Matrixd(3, 4, -1);
    assert(mat1 == [[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]].to!(double[][]));
    
    // duplicating
    auto mdup = matf.dup(); 
    mdup[0][0] = 9;
    assert(mdup == [[9, 2], [3, 4]]);
    assert(matf == [[-15, 2], [3, 4]]);
    
    // diagonal
    auto diag = matDiagonal!float(3, 4, 1);
    assert(diag == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]);
    
    auto diag1 = matDiagonald(3, 4, 1);
    assert(diag1 == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]);
    
    auto diag2 = matDiagonalf(2, 2, 3);
    assert(diag2 == [[3, 0], [0, 3]]);

    // resize
    diag1.resize(2, 3);
    assert(diag1.rows == 2 && diag1.cols == 3);

    diag2.resize(4, 2);
    assert(diag2.rows == 4 && diag2.cols == 2);
    
    diag.resize(2, 2);
    assert(diag == [[1, 0], [0, 0]]);

    // element wise operations
    diag += -1; assert(diag == [[0, -1], [-1, -1]]);
    diag *= 2;  assert(diag == [[0, -2], [-2, -2]]);
    diag /= -2; assert(diag == [[0, 1], [1, 1]]);
    diag -= -1; assert(diag == [[1, 2], [2, 2]]);
    diag--; assert(diag == [[0, 1], [1, 1]]);
    diag++; assert(diag == [[1, 2], [2, 2]]);

    // matrix operations
    diag += diag; assert(diag == [[2, 4], [4, 4]]);
    diag -= diag; assert(diag == [[0, 0], [0, 0]]);
    assert(diag.isZero());
    
    diag++;
    diag = diag + diag; assert(diag == [[2, 2], [2, 2]]);
    diag2.resize(2, 2);
    diag = diag + diag2; assert(diag == [[5, 2], [2, 5]]);
}














