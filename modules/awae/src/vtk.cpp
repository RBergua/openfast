//------------------------------------------------------------------------------
// VTK inflow reading utilities for the AWAE module
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <sstream>

const auto MaxChars{1023};

void convert_string_to_uppercase(std::string &str)
{
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c)
                   { return std::toupper(c); });
}

void copy_string_to_array(const std::string &source, char destination[MaxChars])
{
    const auto n_chars = source.copy(destination, MaxChars);
    for (auto i = n_chars; i < MaxChars; ++i)
    {
        destination[i] = ' ';
    }
}

const auto ErrID_None{0};
const auto ErrID_Info{1};
const auto ErrID_Warn{2};
const auto ErrID_Severe{3};
const auto ErrID_Fatal{4};

extern "C"
{
    void ReadVTK_inflow_info(const char filename[], char desc[MaxChars],
                             int dims[3], double origin[3], double spacing[3],
                             char vec_label[MaxChars], float values[], int *read_values,
                             int *err_stat, char err_msg[MaxChars])
    {
        // Initialize error status and message
        *err_stat = ErrID_Fatal;
        copy_string_to_array("", err_msg);

        // Open the file
        std::ifstream inputFile(filename);

        // If file wasn't opened, return error
        if (!inputFile.is_open())
        {
            copy_string_to_array((std::string("Error opening file: '") + filename + "'"), err_msg);
            return;
        }

        std::string line;
        std::string label;

        std::getline(inputFile, line); // Header
        std::getline(inputFile, line); // Description
        copy_string_to_array(line, desc);

        // Format label
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("ASCII") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find ASCII label", err_msg);
            return;
        }

        // Dataset
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("DATASET") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find DATASET label", err_msg);
            return;
        }
        if (line.find("STRUCTURED_POINTS") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find STRUCTURED_POINTS label", err_msg);
            return;
        }

        // Dimensions
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("DIMENSIONS") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find DIMENSIONS label", err_msg);
            return;
        }
        {
            std::istringstream iss(line);
            iss >> label >> dims[0] >> dims[1] >> dims[2];
        }

        // Origin
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("ORIGIN") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find ORIGIN label", err_msg);
            return;
        }
        {
            std::istringstream iss(line);
            iss >> label >> origin[0] >> origin[1] >> origin[2];
        }

        // Spacing
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("SPACING") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find SPACING label", err_msg);
            return;
        }
        {
            std::istringstream iss(line);
            iss >> label >> spacing[0] >> spacing[1] >> spacing[2];
        }

        // Point data
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("POINT_DATA") == std::string::npos)
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find POINT_DATA label", err_msg);
            return;
        }
        int n_values{0};
        {
            std::istringstream iss(line);
            iss >> label >> n_values;
        }
        if (n_values != (dims[0] * dims[1] * dims[2]))
        {
            copy_string_to_array("Invalid vtk structured_points file: POINT_DATA does not match DIMENSIONS", err_msg);
            return;
        }

        // vector or field data
        std::getline(inputFile, line);
        convert_string_to_uppercase(line);
        if (line.find("VECTORS") != std::string::npos)
        {
            if (line.find("FLOAT") == std::string::npos)
            {
                copy_string_to_array("Invalid VECTORS datatype.  Must be set to float.", err_msg);
                return;
            }
            copy_string_to_array(line.substr(9), vec_label);
        }
        else if (line.find("FIELD") != std::string::npos)
        {
            std::istringstream iss(line);
            int n_arrays{0};
            iss >> label >> label >> n_arrays;
            if (n_arrays != 1)
            {
                copy_string_to_array("Invalid vtk structured_points file: FIELD label must have only 1 array", err_msg);
                return;
            }

            std::getline(inputFile, line);
            convert_string_to_uppercase(line);
            if (line.find("FLOAT") == std::string::npos)
            {
                copy_string_to_array("Invalid FIELD datatype.  Must be set to float.", err_msg);
                return;
            }

            int n_components{0};
            {
                std::istringstream iss(line);
                iss >> label >> n_components >> n_values;
            }
            if (n_components != 3)
            {
                copy_string_to_array("Invalid FIELD components.  Must be set to 3.", err_msg);
                return;
            }
            if (n_values != (dims[0] * dims[1] * dims[2]))
            {
                copy_string_to_array("Invalid vtk structured_points file: FIELD array does not match DIMENSIONS", err_msg);
                return;
            }
        }
        else
        {
            copy_string_to_array("Invalid vtk structured_points file: did not find VECTORS or FIELD label", err_msg);
            return;
        }

        // If reading of values was not requested, return
        if (read_values == 0)
            return;

        // Read values
        for (auto i = 0U; i < n_values; ++i)
        {
            inputFile >> values[i];
        }
    }
}
