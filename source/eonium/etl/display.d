module eonium.etl.display;

/++ DISPLAY module
    - display
+/

import std.conv: to;
import std.stdio: writeln;
import std.string: format;
import std.traits: isNumeric;
import std.array: array;
import std.range: transposed;
import std.algorithm: map, fold, max;

/++
Pretty-prints a 2d array to terminal

Params:
    data = 2d array
    rows = row range to display [from, to] (default: [0, 10])
    cols = column range to display [from, to] (default: [0, max])
    print = print the result to terminal (default: true)

Returns: a formatted string of data in table format

Notes: if rows and cols range are out-of-bounds, then the max array length is used instead for both rows and columns.
+/
string display(T = string)(const T[][] data, const size_t[2] rows = [0, 10], const size_t[2] cols = [0, 0], const bool print = true) 
    in(data !is null && data.length && data[0].length) 
{
    // prepare variables for slicing
    immutable r1 = rows[0];
    immutable r2 = (rows[1] != 0 && rows[1] < data.length ? rows[1] : data.length);
    immutable c1 = cols[0];
    immutable c2 = (cols[1] != 0 && cols[1] < data[0].length ? cols[1] : data[0].length);
    
    // padding
    immutable totalRowPadding = data.length.to!string.length;
    immutable totalColPadding = data[0].length.to!string.length; 
    

    // create a copy we can freely modify
    const df = data[r1..r2].dup()       // duplicate the array within the rows range
        .to!(string[][])                // convert it to string
        .map!(a => a[c1..c2])           // extract columns within the cols range
        .array;                         // convert from the map's result type back to df array type
    
    // calculates the maximum element-wise padding of the array (the length of the widest array element)
    string calcPadding(S)(S arr) {
        return arr.map!(a => a.length + totalColPadding + totalRowPadding)  // get elements length (adjusted with padding)
            .fold!max                   // find max length
            .to!string;                 // convert it to string
    }
    
    // max elemnt-wise padding for each column
    const padding = df.dup().transposed().map!(a => calcPadding(a)).array;
    
    // formatted string
    string fmt = format("\nRow/Column\nIndex\n%" ~ totalRowPadding.to!string  ~ "s%4s", "-", " ");

    // format column index
    foreach(i; c1..c2) {
        fmt ~= format("%" ~ padding[i] ~ "s%4s", i, " ");
    }
    fmt ~= "\n";

    // format row index and the data itself
    foreach(i, row; df) {
        fmt ~= format("%" ~ totalRowPadding.to!string ~ "s%4s", i + r1, " ");
        foreach(j, entry; row) {
            fmt ~= format("%" ~ padding[j] ~ "s%4s", entry, " ");
        }
        fmt ~= "\n";
    }
    
    // format basic information about the data
    fmt ~= "\n"
        ~ format("dim(displayed): [ %" ~ totalRowPadding.to!string ~ "s x %" ~ totalColPadding.to!string ~ "s ]\n", r2 - r1, c2 - c1)
        ~ format("dim(dataframe): [ %" ~ totalRowPadding.to!string ~ "s x %" ~ totalColPadding.to!string ~ "s ]\n", data.length, data[0].length)
        ~ "\n";

    // print the formatted string
    if(print) {
        fmt.writeln; 
    }

    return fmt;
}

/* --------------------- UNITTESTS --------------------- */

unittest {
    import eonium.etl.csv;
    import eonium.la.mat;
    
    // display array
    float[][] arr = [
        [1, 21], 
        [3, 4], 
        [155, 0]
    ];
    //arr.display;

    // display from csv #1
    auto csvData = "data/csv/test_csvRead2.csv".csvRead;
    //csvData.display;
    
    // display from matrix
    auto mat = Matrixf(arr);
    //mat.display;

    // display from csv #2
    auto csvData2 = "data/csv/test_stud.csv".csvRead;
    auto str = csvData2.display([0, 0], [0, 6], false); // display all rows and first 6 cols
    //str.writeln;

    //auto d = "data/csv/2016.csv".csvRead(",");
    //d.display([0, 0], [0, 0]); // display all rows and cols
}













