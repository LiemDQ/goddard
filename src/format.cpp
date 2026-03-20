#include "goddard/format.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace Goddard {

static std::string right_align(const std::string& s, int width) {
    if (static_cast<int>(s.size()) >= width) return s;
    return std::string(width - s.size(), ' ') + s;
}

static std::string left_align(const std::string& s, int width) {
    if (static_cast<int>(s.size()) >= width) return s;
    return s + std::string(width - s.size(), ' ');
}

std::string format_fixed(double value, int width, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return right_align(oss.str(), width);
}

std::string format_cea_engineering(double value, int width) {
    if (value == 0.0) {
        return right_align("0.0000 0", width);
    }

    bool negative = value < 0.0;
    double abs_val = std::abs(value);

    int exponent = static_cast<int>(std::floor(std::log10(abs_val)));
    double mantissa = abs_val / std::pow(10.0, exponent);

    // Handle rounding edge case: mantissa could round up to 10.0
    if (mantissa >= 9.99995) {
        mantissa /= 10.0;
        exponent += 1;
    }

    std::ostringstream oss;
    if (negative) oss << "-";
    oss << std::fixed << std::setprecision(4) << mantissa;
    if (exponent >= 0) {
        oss << " " << exponent;
    } else {
        oss << exponent;
    }
    return right_align(oss.str(), width);
}

std::string format_mass_fraction(double value, int width) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(5) << value;
    return right_align(oss.str(), width);
}

TextTable::TextTable(int label_width, int column_width)
    : m_label_width(label_width), m_column_width(column_width) {}

void TextTable::set_headers(const std::vector<std::string>& headers) {
    m_rows.push_back(Row{RowType::HEADER, "", headers});
}

void TextTable::add_row(const std::string& label, const std::vector<std::string>& values) {
    m_rows.push_back(Row{RowType::DATA, label, values});
}

void TextTable::add_blank_line() {
    m_rows.push_back(Row{RowType::BLANK, "", {}});
}

void TextTable::add_section_header(const std::string& header) {
    m_rows.push_back(Row{RowType::SECTION_HEADER, header, {}});
}

std::string TextTable::render() const {
    std::string result;

    for (const auto& row : m_rows) {
        switch (row.type) {
            case RowType::HEADER: {
                result += std::string(m_label_width, ' ');
                for (const auto& h : row.values) {
                    result += right_align(h, m_column_width);
                }
                result += "\n";
                break;
            }
            case RowType::DATA: {
                result += left_align(row.label, m_label_width);
                for (const auto& v : row.values) {
                    result += right_align(v, m_column_width);
                }
                result += "\n";
                break;
            }
            case RowType::BLANK: {
                result += "\n";
                break;
            }
            case RowType::SECTION_HEADER: {
                result += row.label + "\n";
                break;
            }
        }
    }

    return result;
}

} // namespace Goddard
