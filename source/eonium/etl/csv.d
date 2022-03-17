module eonium.etl.csv;

/++ CSV module
    - csvRead
    - csvWrite
+/

import std.stdio: File, writeln, stderr;
import std.file: exists;
import std.conv: to;
import std.traits: isNumeric;
import std.string: format;
import std.array: empty, split, array, replace;
import std.algorithm: map, each, joiner, remove, among;

/++
Read a CSV file into memory
Params:
    filename = path to the file
    sep = seperator, (default: ';')
    header = does CSV contain a header (default: true)
    preallocate = pre-allocate N number of rows (default: 100)
Returns: 'string[][]' upon success, 'null' upon failure
Note: if a CSV file contains less/more entries than in the first row, empty entries are appended/removed.
+/
T[][] csvRead(T = string)(const string filename, const string sep = ";", bool header = true, const size_t preallocate = 100) in {
    assert(filename.exists, format("<%s> does not exist!", filename));
} do {
    // opening a file
    File file = File(filename, "r");
    scope(exit) { file.close(); }

    // check if file was opened
    if(!file.isOpen) {
        stderr.writeln(format("Cannot open <%s>!", filename));
        return null;
    }

    // reading data from the file
    size_t previousEntriesLen = 0;
    T[][] data; data.reserve(preallocate);
    foreach(record; file.byLine.map!(row => row.replace("\"", "").split(sep)).array.remove!(row => row.empty)) {
        // if header is present, skip
        if(header) {
            header = false;
            continue;
        }

        // save number of entries in the first row
        if(previousEntriesLen == 0) {
            previousEntriesLen = record.length;
        }

        // if CSV file is damaged, try to fix it
        if(record.length < previousEntriesLen) {
            // add empty entries
            while(record.length < previousEntriesLen) {
                record ~= "".to!(char[]);
            }
        } else {
            // remove entries
            record = record[0..previousEntriesLen];
        }

        // save the row
        try {
            data ~= record.map!(entry => (entry.empty || entry.among("N/A", "n/a", "na")) ? (isNumeric!T ? "nan" : "NA") : (entry)).array.to!(T[]);
        } catch(Exception e) {
            stderr.writeln(format("<%s> file is damaged. Unable to repair the CSV file. Error: %s", filename, e.msg));
        }
    }

    return data;
}

/++
Writes data to a CSV file
Params:
    data = 2D array of type T, (default: string[][])
    filename = file name
    header = header names seperated by 'sep' delimiter (default: null)
    sep = seperator (default: ';')
    quoted = quote each entry with double quotes (default: true)
Note: the header separator must match 'sep'.
+/
void csvWrite(T = string)(const T[][] data, const string filename, const string header = null, const string sep = ";", const bool quoted = true) {
    // if data is empty, do nothing
    if(data.empty) {
        stderr.writeln("Data is empty, nothing to write!");
        return;
    }

    // open file
    File file = File(filename, "w");
    scope(exit) { file.close(); }
    if(!file.isOpen) {
        stderr.writeln(format("Could not open <%s>!", filename));
        return;
    }

    // quote data
    immutable quote = quoted ? "\"" : "";

    // write header
    if(!header.empty) {
        file.writeln(header.split(sep).map!(entry => format("%s%s%s", quote, entry, quote)).joiner(sep));
    }

    // write csv contents
    data.each!(row => file.writeln(row.map!(entry => format("%s%s%s", quote, entry, quote)).joiner(sep)));
}

/* --------------------- UNITTESTS --------------------- */

unittest {
    import std;

    auto data = "data/csv/test_csvRead.csv".csvRead;
    assert(data.length == 3);

    auto data2 = "data/csv/test_csvRead2.csv".csvRead;
    assert(data2.length == 7);

    data2.csvWrite("data/csv/test_csvWrite.csv");
    assert("data/csv/test_csvWrite.csv".exists);
}