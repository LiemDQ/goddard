#pragma once

#include <string>
#include <vector>

namespace Goddard {

/// Format a number in fixed-point notation, right-aligned to the given width.
std::string format_fixed(double value, int width, int precision);

/// Format a number in CEA engineering notation.
/// CEA uses: "9.0033 0" for 9.0033e0, "1.4339-1" for 1.4339e-1.
/// Mantissa has 4 decimal places; exponent is space+digit (positive) or -digit (negative).
std::string format_cea_engineering(double value, int width);

/// Format a mass fraction value (5 decimal places, right-aligned).
std::string format_mass_fraction(double value, int width);

/// A simple column-aligned text table for CEA-style report output.
/// Values are pre-formatted as strings before being added.
/// The table handles only layout and alignment.
class TextTable {
public:
    TextTable(int label_width = 18, int column_width = 10);

    /// Set column headers (e.g. {"CHAMBER", "THROAT", "EXIT", "EXIT"}).
    void set_headers(const std::vector<std::string>& headers);

    /// Add a data row with a label and pre-formatted values.
    void add_row(const std::string& label, const std::vector<std::string>& values);

    /// Add a blank line.
    void add_blank_line();

    /// Add a section header line (e.g. "PERFORMANCE PARAMETERS").
    void add_section_header(const std::string& header);

    /// Render the complete table to a string.
    std::string render() const;

private:
    int m_label_width;
    int m_column_width;

    enum class RowType { HEADER, DATA, BLANK, SECTION_HEADER };

    struct Row {
        RowType type;
        std::string label;
        std::vector<std::string> values;
    };

    std::vector<Row> m_rows;
};

} // namespace Goddard
