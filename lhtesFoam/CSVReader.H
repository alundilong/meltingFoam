#ifndef CSVREADER_H
#define CSVREADER_H

#include <string>
#include <vector>
#include <iostream>

class CSVReader {
public:
    // Structure to hold a row of data
    struct DataRow {
        double scalar;
        int vtkOriginalPointId;
        double point0, point1, point2;
    };

    // Constructor
    explicit CSVReader(const std::string& filename);

    // Method to load the CSV file
    bool loadCSV();

    // Method to print the data
    void printData() const;

    // Getter for data rows
    const std::vector<DataRow>& getData() const;

private:
    std::string filename;
    std::vector<DataRow> dataRows;
};

#endif // CSVREADER_H

