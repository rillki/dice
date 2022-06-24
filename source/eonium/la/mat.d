module eonium.la.mat;

/++ MATRIX module
    - struct Matrix
        - Matrix: from 2d array
        - Matrix: from size (rows, cols)
        - dup
        - resize
        - transpose
        - inverse
    - Identity
+/

import std.math: abs, sqrt, mround = round;
import std.traits: isFloatingPoint;
import std.algorithm.searching: canFind;

alias Matrixf = Matrix!float;
alias Matrixd = Matrix!double;
alias Identityf = Identity!float;
alias Identityd = Identity!double;

struct Matrix(T = float) if(isFloatingPoint!T) {
    private {
        T[] data;
        size_t r;
        size_t c;
    }
    
    T[][] m;
    alias m this;
    
    /// Construct matrix from 2d array
    this(const T[][] arr) in (arr !is null && arr.length && arr[0].length) {
        // save matrix size
        r = arr.length;
        c = arr[0].length;
    
        // copy matrix data
        data.length = r * c;
        foreach(i, row; arr) {
            foreach(j, val; row) {
                data[i * c + j] = val;
            }
        }

        // set slices
        updateSlices();
    }

    /// Construct matrix from 1d array
    this(const T[] arr, const size_t r, const size_t c) in(arr !is null && arr.length && r && c && arr.length <= r * c) {
        // save matrix size
        this.r = r;
        this.c = c;

        // copy values
        data.length = this.r * this.c;
        data[] = arr[];

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
    
    /// Checks if matrix is a square matrix
    bool isSquare() const {
        return (r == c);
    }

    /// Returns a copy of a matrix
    Matrix!T dup() const {
        return Matrix!T(this);
    }

    /++ 
    Resize matrix
    
    Params:
        r = rows
        c = columns
        stable = keep values position (default: false)
    +/
    void resize(const size_t r, const size_t c, const bool stable = false) in (r && c) {
        if(stable) {
            auto mcopy = Matrix!T(r, c);
            foreach(i; 0..mcopy.rows) {
                foreach(j; 0..mcopy.cols) {
                    mcopy[i][j] = this[i][j];
                }
            }

            this = mcopy;
        } else {
            this.r = r;
            this.c = c;
            
            data.length = this.r * this.c;
            updateSlices();
        }
    }

    /++
    Drops matrix rows and columns
    
    Params:
        indexes = index number
        axis = { 0: rows, 1: cols } (default: 0)
    +/
    void drop(const size_t[] indexes, const bool axis = 0) in (indexes !is null) {
        Matrix!T mat;
        
        // dropping rows
        if(axis == 0) {
            // create a new empty matrix
            mat = Matrix!T(this.rows - indexes.length, this.cols);

            size_t offset = 0;
            foreach(i; 0..this.rows) {
                // drop rows
                if(indexes.canFind(i)) { offset++; continue; }
                
                // copy rows that we don't drop
                mat[i - offset] = this[i];
            }
        } else { // droping rows
            // create a new empty matrix
            mat = Matrix!T(this.rows, this.cols - indexes.length);
            
            foreach(i; 0..this.rows) {
                size_t offset = 0;
                foreach(j; 0..this.cols) {
                    // drop columns
                    if(indexes.canFind(j)) { offset++; continue; }
                    // copy column values that we don't drop
                    mat[i][j - offset] = this[i][j];
                }
            }
        }
        
        // save it
        this = mat;
    }

    /+ ---------- ELEMENT-WISE OPERATIONS ---------- +/

    /// Element wise increment/decrement (++, --)
    ref Matrix!T opUnary(string op)() if(op == "++" || op == "--") {
        mixin(op ~ "data[];");
        return this;
    }
    
    /// Element wise value operations (+=, -=, *=, /=)
    ref Matrix!T opOpAssign(string op)(const T value) if(op == "+" || op == "-" || op == "*" || op == "/") {
        mixin("data[] " ~ op ~ "= value;");
        return this;
    }

    /// Matrix element wise operations (+, -, *)
    Matrix!T opBinary(string op)(const T value) const if(op == "+" || op == "-" || op == "*") {
        return this.dup().opOpAssign!op(value);
    }

    /// Matrix element wise operations (+, -, *)
    Matrix!T opBinaryRight(string op)(const T value) const if(op == "+" || op == "-" || op == "*") {
        return this.dup().opOpAssign!op(value);
    }

    /+ ---------- MATRIX-WISE OPERATIONS ---------- +/
    
    /// Matrix-wise operations (+=, -=)
    ref Matrix!T opOpAssign(string op)(const Matrix!T mat) if(op == "+" || op == "-") in(this.r == mat.rows && this.c == mat.cols) {
        mixin("data[] " ~ op ~ "= mat.array[];");
        return this;
    } 
    
    /// Matrix-wise operations (*=)
    ref Matrix!T opOpAssign(string op)(const Matrix!T mat) if(op == "*") in(this.r == mat.cols) {
        auto ret = Matrix!T(this.r, mat.cols, 0);
        
        // matrix multiplication
        foreach(i; 0..r) {
            foreach(j; 0..mat.cols) {
                foreach(k; 0..mat.rows) {
                    ret[i][j] += this[i][k] * mat[k][j];
                }
            }
        }
        
        // update matrix
        this = ret;

        return this;
    }

    /// Matrix-wise operations (+, -, *)
    Matrix!T opBinary(string op)(const Matrix!T mat) const if(op == "+" || op == "-" || op == "*") {
        return this.dup().opOpAssign!op(mat);
    }
    
    /// Matrix inplace transposition
    ref Matrix!T transpose() {
        // transpose the matrix (from RosettaCode, the C version, uses permutations)
        // RosettaCode: http://www.rosettacode.org/wiki/Matrix_transposition#C
        // Wiki article on permutations: https://en.wikipedia.org/wiki/In-place_matrix_transposition
        for(size_t start = 0; start < r * c; start++) {
            size_t next = start;
            size_t i = 0;

            do {
                i++;
                next = (next % r) * c + next / r;
            } while (next > start);

            if(next < start || i == 1) {
                continue;
            }

            immutable T tmp = data[next = start];
            do {
                i = (next % r) * c + next / r;
                data[next] = (i == start) ? tmp : data[i];
                next = i;
            } while (next > start);
        }

        // update matrix size
        immutable tmp = r;
        r = c;
        c = tmp;

        // update slices pointers
        updateSlices();

        return this;
    }

    /// Inplace matrix inverse
    /// Notes: doesn't do anything if matrix is singular (from RosettaCode, the C version, using Gauss-Jordan algorithm)
    ref Matrix!T inverse() in(r == c) {
        // RosettaCode: http://www.rosettacode.org/wiki/Gauss-Jordan_matrix_inversion#C
        auto mat = Identity!T(r, c, 1);

        // Frobenius norm of a
        T f = 0;
        T g = 0;
        for(size_t i = 0; i < r; ++i) {
            for(size_t j = 0; j < r; ++j) {
                g = data[j + i * r];
                f += g * g;
            }
        }

        f = sqrt(f);
        T tol = f * 2.2204460492503131e-016;

        // main loop
        for (size_t k = 0; k < r; ++k) {
            // find pivot
            size_t p = k;
            f = abs(data[k + k * r]);
            for (size_t i = k + 1; i < r; ++i) {
                g = abs(data[k + i * r]);
                if (g > f) {
                    f = g;
                    p = i;
                }
            }

            // matrix is singluar
            if (f < tol) {
                return this;
            }

            // swap rows
            if (p != k) {
                for (size_t j = k; j < r; ++j) {
                    f = data[j + k * r];
                    data[j + k * r] = data[j + p * r];
                    data[j + p * r] = f;
                }

                for (size_t j = 0; j < r; ++j) {
                    f = mat[k][j];
                    mat[k][j] = mat[p][j];
                    mat[p][j] = f;
                }
            }

            // scale row so pivot is 1
            f = 1.0 / data[k + k * r];
            for (size_t j = k; j < r; ++j) {
                data[j + k * r] *= f;
            }

            for (size_t j = 0; j < r; ++j) {
                mat[k][j] *= f;
            }

            // Subtract to get zeros
            for (size_t i = 0; i < r; ++i) {
                if (i == k) {
                    continue;
                }

                f = data[k + i * r];
                for (size_t j = k; j < r; ++j) {
                    data[j + i * r] -= data[j + k * r] * f;
                }

                for (size_t j = 0; j < r; ++j) {
                    mat[i][j] -= mat[k][j] * f;
                }
            }
        }

        // update matrix
        this = mat;

        return this;
    }

    /++ 
    Inplace round matrix elements
    
    Params:
        n = number of digits after the floating point '.'
    
    Returns: a rounded Matrix!T
    +/
    ref Matrix!T round(const size_t n = 2) {
        immutable scaler = 10^^n;
        foreach(ref d; data) {
            d *= scaler;
            d = mround(d);
            d /= scaler;
        }

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
Creates an identity matrix with custom value filled along the diagonal (default: 1)

Params:
    rows = number of rows
    cols = number of cols
    value = diagonal fill value

Returns: Matrix!T
+/
Matrix!T Identity(T = float)(const size_t r, const size_t c, const T value = 1) if(isFloatingPoint!T) in (r && c) {
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
    auto diag = Identity!float(3, 4, 1);
    assert(diag == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]);
    
    auto diag1 = Identityd(3, 4, 1);
    assert(diag1 == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]);
    
    auto diag2 = Identityf(2, 2, 3);
    assert(diag2 == [[3, 0], [0, 3]]);

    // resize
    diag1.resize(2, 3);
    assert(diag1.rows == 2 && diag1.cols == 3);

    diag2.resize(4, 2);
    assert(diag2.rows == 4 && diag2.cols == 2);
    
    diag.resize(2, 2);
    assert(diag == [[1, 0], [0, 0]]);

    // ---- ELEMENT-WISE OPERATIONS ----

    // +=, -=, *=, /=, ++, --
    diag += -1; assert(diag == [[0, -1], [-1, -1]]);
    diag *= 2;  assert(diag == [[0, -2], [-2, -2]]);
    diag /= -2; assert(diag == [[0, 1], [1, 1]]);
    diag -= -1; assert(diag == [[1, 2], [2, 2]]);
    diag--; assert(diag == [[0, 1], [1, 1]]);
    diag++; assert(diag == [[1, 2], [2, 2]]);
    --diag; assert(diag == [[0, 1], [1, 1]]);
    ++diag; assert(diag == [[1, 2], [2, 2]]);

    // +, -, *
    auto diagCopy = diag.dup();
    diag = diagCopy + 2;   assert(diag == [[3, 4], [4, 4]]);
    diag = diag - 1;       assert(diag == [[2, 3], [3, 3]]);
    diag = diagCopy * 0.5; assert(diag == [[0.5, 1.0], [1.0, 1.0]]);
    diag = 3 * diag;       assert(diag == [[1.5, 3.0], [3.0, 3.0]]);
    
    // ---- MATRIX-WISE OPERATIONS ---

    // +=, -=
    diag += diag; assert(diag == [[3, 6], [6, 6]]);
    diag -= diag; assert(diag.isZero());
    
    // *=
    diag += diagCopy; assert(diag == [[1, 2], [2, 2]]);
    diag *= diag;     assert(diag == [[5, 6], [6, 8]]);
    
    // +, -
    diag = diag + diagCopy; assert(diag == [[6, 8], [8, 10]]);
    diag = diagCopy - diag; assert(diag == [[-5, -6], [-6, -8]]);
    
    // *
    auto m1 = Matrixf([[2, 3], [-5, 6], [9, -7]]);
    auto m2 = Matrixf([[1, -2, 0], [3, 4, -5]]);
    auto m3 = m1 * m2; assert(m3 == [[11, 8, -15], [13, 34, -30], [-12, -46, 35]]);
    
    // matrix transposition
    auto m1T = m1.dup().transpose(); assert(m1T == [[2, -5, 9], [3, 6, -7]]);
    m2.transpose; assert(m2 == [[1, 3], [-2, 4], [0, -5]]);

    // inverse
    auto m4 = Matrixf([
        [-1, -2, 3, 2],
        [-4, -1, 6, 2],
        [7, -8, 9, 1],
        [1, -2, 1, 3]
    ]);
    
    auto m4Inv = m4.dup().inverse(); assert(m4Inv != m4);
    assert(m4.inverse.inverse == m4);
    
    // identity inverse
    auto I = Identityf(3, 3);
    assert(I.inverse == I);

    // resize(stable = true)
    auto m4dup = m4.dup();

    m4.resize(3, 2);
    assert(m4.rows == 3 && m4.cols == 2);

    m4dup.resize(3, 2, true);
    assert(m4dup.rows == 3 && m4dup.cols == 2);
    assert(m4dup.round() == [[-1, -2], [-4, -1], [7, -8]]);

    // drop rows
    auto m5 = Matrixf([
        [-1, -2, 3, 2],
        [-4, -1, 6, 2],
        [7, -8, 9, 1],
        [1, -2, 1, 3]
    ]);
    
    auto mdropRow = m5.dup; mdropRow.drop([0, 2]);
    assert(mdropRow == [[-4, -1, 6, 2], [1, -2, 1, 3]]);

    // ensure that m5 value does not change when mdrop is modified
    mdropRow[0][0] = 153;
    assert(m5[1][0] == -4);

    // drop columns
    auto mdropCol = mdropRow.dup; mdropCol.drop([2, 3], 1);
    assert(mdropCol == [[153, -1], [1, -2]].to!(float[][]));
}














