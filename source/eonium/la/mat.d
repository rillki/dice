module eonium.la.mat;

/++ MATRIX module
    - struct Matrix
        - Matrix: from 2d array
        - Matrix: from size (rows, cols)
        - dup
    - matDiagonal
+/

import std.algorithm: copy;
import std.traits: isFloatingPoint;

alias Matrixf = Matrix!float;
alias Matrixd = Matrix!double;
alias matDiagonalf = matDiagonal!float;
alias matDiagonald = matDiagonal!double;

struct Matrix(T = float) if(isFloatingPoint!T) {
    T[][] m;
    T[] data;
    size_t rows;
    size_t cols;
    alias m this;
    
    /// Construct matrix from 2d array
    this(const T[][] arr) in (arr.length && arr[0].length) {
        // save matrix size
        rows = arr.length;
        cols = arr[0].length;

        // copy matrix data
        data.length = rows * cols;
        foreach(i, row; arr) {
            foreach(j, val; row) {
                data[i * rows + j] = val;
            }
        }

        // set slices
        updateSlices();
    }

    /// Construct a matrix from size (rows, cols) and fill with value
    this(const size_t r, const size_t c, const T value = 0) in (r && c) {
        // save matrix size
        rows = r;
        cols = c;

        // create an empty matrix and fill it with value
        data.length = rows * cols;
        data[] = value;
        
        // set slices
        updateSlices();
    }
    
    /// Returns a copy of a matrix
    Matrix!T dup() const {
        return Matrix!T(this);
    }
    
    private {
        /// Updates slice pointers
        void updateSlices() {
            m.length = rows;
            foreach(i, ref row; m) {
                immutable s = i*cols;
                row = data[s..(s + cols)];
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
Matrix!T matDiagonal(T = float)(const size_t rows, const size_t cols, const T value = 1) if(isFloatingPoint!T) {
    auto mat = Matrix!T(rows, cols, 0);
    for(size_t i = 0; i < (rows < cols ? rows : cols); i++) {
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
    
    auto diag1 = matDiagonalf(3, 4, 1);
    assert(diag1 == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]);
    
    auto diag2 = matDiagonald(2, 2, 3);
    assert(diag2 == [[3, 0], [0, 3]]);
}















