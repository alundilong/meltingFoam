#include "CSVReader.H"
#include <fstream>
#include <sstream>
#include <string>

// Constructor
CSVReader::CSVReader(const std::string& filename) : filename(filename) {}

// Method to load the CSV file
bool CSVReader::loadCSV() {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;

    // Skip the header line
    if (!std::getline(file, line)) {
        std::cerr << "Error: File is empty or header missing." << std::endl;
        return false;
    }

    // Read each line of the CSV file
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        DataRow row;

        // Parse each value
        std::getline(ss, token, ','); // Scalars
        row.scalar = std::stod(token);

        std::getline(ss, token, ','); // vtkOriginalPointIds
        row.vtkOriginalPointId = std::stoi(token);

        std::getline(ss, token, ','); // Points:0
        row.point0 = std::stod(token);

        std::getline(ss, token, ','); // Points:1
        row.point1 = std::stod(token);

        std::getline(ss, token, ','); // Points:2
        row.point2 = std::stod(token);

        // Store the parsed data
        dataRows.push_back(row);
    }

    file.close();
    return true;
}

// Method to print the data
void CSVReader::printData() const {
    for (const auto& row : dataRows) {
        std::cout << "Scalar: " << row.scalar
                  << ", vtkOriginalPointId: " << row.vtkOriginalPointId
                  << ", Points: (" << row.point0 << ", " << row.point1 << ", " << row.point2 << ")"
                  << std::endl;
    }
}

// Getter for data rows
const std::vector<CSVReader::DataRow>& CSVReader::getData() const {
    return dataRows;
}

