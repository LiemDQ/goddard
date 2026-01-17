#!/usr/bin/env python3
"""
CEA Parser - Extract test data from NASA CEA output files

This tool parses CEA output files and extracts key thermodynamic properties
into structured JSON format for use in Goddard unit tests.
"""

import re
import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, asdict


@dataclass
class CEAConditions:
    """Input conditions from CEA file"""
    pressure_psia: float
    of_ratio: float
    fuel: str
    oxidizer: str
    fuel_temp: float
    fuel_energy: float  # KJ/KG-MOL
    oxidizer_temp: float
    oxidizer_energy: float  # KJ/KG-MOL
    area_ratios: List[float]  # Supersonic area ratios
    case_name: str
    problem_type: str


@dataclass
class CEAState:
    """Thermodynamic state at a specific condition"""
    location: str  # CHAMBER, THROAT, EXIT
    pressure_ratio: float
    pressure_bar: float
    temperature_k: float
    density_kg_m3: float
    enthalpy_kj_kg: float
    internal_energy_kj_kg: float
    gibbs_kj_kg: float
    entropy_kj_kg_k: float
    molecular_weight: float
    dlv_dlp_t: float
    dlv_dlt_p: float
    cp_kj_kg_k: float
    gamma: float
    sound_speed_m_s: float
    mach_number: float
    mass_fractions: Dict[str, float]


@dataclass
class CEAPerformance:
    """Rocket performance parameters"""
    area_ratio: float  # Ae/At
    cstar_m_s: float
    cf: float
    isp_m_s: float
    ivac_m_s: float


@dataclass
class CEAResult:
    """Complete CEA analysis result"""
    conditions: CEAConditions
    equilibrium_states: List[CEAState]
    frozen_states: List[CEAState]
    equilibrium_performance: List[CEAPerformance]
    frozen_performance: List[CEAPerformance]


class CEAParser:
    """Parser for NASA CEA output files"""
    
    def __init__(self):
        self.reactant_pattern = re.compile(
            r'(FUEL|OXIDANT)\s+(\w+(?:\([LG]\))?)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)'
        )
        self.of_pattern = re.compile(r'O/F=\s*(\d+\.\d+)')
        self.pressure_pattern = re.compile(r'p,psia=\s*(\d+(?:\.\d+)?)')
        
    def parse_file(self, filepath: Path) -> CEAResult:
        """Parse a CEA output file and return structured data"""
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Try to find corresponding input file for more complete data
        input_filepath = self._find_input_file(filepath)
        input_conditions = None
        if input_filepath and input_filepath.exists():
            input_conditions = self._parse_input_file(input_filepath)
        
        # Parse conditions from output file
        output_conditions = self._parse_conditions(content)
        
        # Merge input and output conditions (input takes precedence for missing data)
        conditions = self._merge_conditions(input_conditions, output_conditions)
        
        # Split equilibrium and frozen sections
        eq_section, frozen_section = self._split_sections(content)
        
        # Parse equilibrium results
        eq_states = self._parse_states(eq_section, "equilibrium")
        eq_performance = self._parse_performance(eq_section)
        # Parse frozen results
        frozen_states = self._parse_states(frozen_section, "frozen")
        frozen_performance = self._parse_performance(frozen_section)
        
        return CEAResult(
            conditions=conditions,
            equilibrium_states=eq_states,
            frozen_states=frozen_states,
            equilibrium_performance=eq_performance,
            frozen_performance=frozen_performance
        )
    
    def _parse_conditions(self, content: str) -> CEAConditions:
        """Extract input conditions from CEA output"""
        # Find pressure
        pressure_match = self.pressure_pattern.search(content)
        pressure = float(pressure_match.group(1)) if pressure_match else 0.0
        
        # Find O/F ratio
        of_match = self.of_pattern.search(content)
        of_ratio = float(of_match.group(1)) if of_match else 0.0
        
        # Find reactants
        reactants = self.reactant_pattern.findall(content)
        
        fuel_name = fuel_temp = fuel_energy = ""
        oxidizer_name = oxidizer_temp = oxidizer_energy = ""
        
        for reactant_type, name, wt_frac, energy, temp in reactants:
            if reactant_type == "FUEL":
                fuel_name = name
                fuel_temp = float(temp)
                fuel_energy = float(energy)
            elif reactant_type == "OXIDANT":
                oxidizer_name = name
                oxidizer_temp = float(temp)
                oxidizer_energy = float(energy)
        
        return CEAConditions(
            pressure_psia=pressure,
            of_ratio=of_ratio,
            fuel=fuel_name,
            oxidizer=oxidizer_name,
            fuel_temp=fuel_temp,
            fuel_energy=fuel_energy,
            oxidizer_temp=oxidizer_temp,
            oxidizer_energy=oxidizer_energy,
            area_ratios=[],
            case_name="",
            problem_type=""
        )
    
    def _find_input_file(self, output_filepath: Path) -> Optional[Path]:
        """Find corresponding input file for an output file"""
        base_name = output_filepath.stem
        input_path = output_filepath.parent / f"{base_name}.input"
        return input_path if input_path.exists() else None
    
    def _parse_input_file(self, input_filepath: Path) -> CEAConditions:
        """Parse CEA input file to extract complete problem setup"""
        with open(input_filepath, 'r') as f:
            content = f.read()
        
        # Extract basic info
        pressure_match = self.pressure_pattern.search(content)
        pressure = float(pressure_match.group(1)) if pressure_match else 0.0
        
        of_match = self.of_pattern.search(content)
        of_ratio = float(of_match.group(1)) if of_match else 0.0
        
        # Extract case name
        case_match = re.search(r'case=(\w+)', content)
        case_name = case_match.group(1) if case_match else ""
        
        # Extract problem type
        prob_match = re.search(r'# Problem Type: "([^"]*)"', content)
        problem_type = prob_match.group(1) if prob_match else ""
        
        # Extract area ratios
        area_ratios = []
        supar_match = re.search(r'supar=\s*([\d\s,\.]+)', content)
        if supar_match:
            area_str = supar_match.group(1)
            for part in area_str.replace(',', ' ').split():
                try:
                    area_ratios.append(float(part))
                except ValueError:
                    continue
        
        # Extract reactants
        fuel_name = fuel_temp = fuel_energy = ""
        oxidizer_name = oxidizer_temp = oxidizer_energy = ""
        
        # Look for fuel line
        fuel_match = re.search(r'fuel\s+(\S+).*?t,k=\s*(\d+\.\d+)', content)
        if fuel_match:
            fuel_name = fuel_match.group(1)
            fuel_temp = float(fuel_match.group(2))
        else:
            # Try simpler pattern
            fuel_match = re.search(r'fuel\s+(\S+)', content)
            if fuel_match:
                fuel_name = fuel_match.group(1)
                fuel_temp = 298.15  # Default
        
        # Look for oxidizer line
        oxid_match = re.search(r'oxid\s+(\S+).*?t,k=\s*(\d+\.\d+)', content)
        if oxid_match:
            oxidizer_name = oxid_match.group(1)
            oxidizer_temp = float(oxid_match.group(2))
        else:
            # Try simpler pattern
            oxid_match = re.search(r'oxid\s+(\S+)', content)
            if oxid_match:
                oxidizer_name = oxid_match.group(1)
                oxidizer_temp = 298.15  # Default
        
        return CEAConditions(
            pressure_psia=pressure,
            of_ratio=of_ratio,
            fuel=fuel_name,
            oxidizer=oxidizer_name,
            fuel_temp=fuel_temp,
            fuel_energy=fuel_energy,
            oxidizer_temp=oxidizer_temp,
            oxidizer_energy=oxidizer_energy,
            area_ratios=area_ratios,
            case_name=case_name,
            problem_type=problem_type
        )
    
    def _merge_conditions(self, input_conditions: Optional[CEAConditions], 
                         output_conditions: CEAConditions) -> CEAConditions:
        """Merge input and output conditions, preferring input data when available"""
        if not input_conditions:
            return output_conditions
        
        # Use input conditions as base, fill in missing data from output
        return CEAConditions(
            pressure_psia=input_conditions.pressure_psia or output_conditions.pressure_psia,
            of_ratio=input_conditions.of_ratio or output_conditions.of_ratio,
            fuel=input_conditions.fuel or output_conditions.fuel,
            oxidizer=input_conditions.oxidizer or output_conditions.oxidizer,
            fuel_temp=input_conditions.fuel_temp or output_conditions.fuel_temp,
            fuel_energy=output_conditions.fuel_energy,  # Usually only in output
            oxidizer_temp=input_conditions.oxidizer_temp or output_conditions.oxidizer_temp,
            oxidizer_energy=output_conditions.oxidizer_energy,  # Usually only in output
            area_ratios=input_conditions.area_ratios,
            case_name=input_conditions.case_name,
            problem_type=input_conditions.problem_type
        )
    
    def _split_sections(self, content: str) -> tuple[str, str]:
        """Split content into equilibrium and frozen sections"""
        frozen_start = content.find("THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN")
        if frozen_start == -1:
            return content, ""
        
        equilibrium_section = content[:frozen_start]
        frozen_section = content[frozen_start:]
        
        return equilibrium_section, frozen_section
    
    def _parse_states(self, section: str, section_type: str) -> List[CEAState]:
        """Parse thermodynamic states from a section"""
        states = []
        lines = section.split('\n')
        
        # Find header line with CHAMBER THROAT EXIT
        header_idx = -1
        for i, line in enumerate(lines):
            if "CHAMBER" in line and "THROAT" in line and "EXIT" in line:
                header_idx = i
                break
        
        if header_idx == -1:
            return states
        
        header_line = lines[header_idx].strip()
        locations = header_line.split()  # Standard CEA format
        
        
        # Parse data table
        data_dict = {}
        param_mapping = {
            'Pinf/P': 'pinf_p_ratio',
            'P, BAR': 'pressure_bar', 
            'T, K': 'temperature_k',
            'RHO, KG/CU M': 'density_kg_m3',
            'H, KJ/KG': 'enthalpy_kj_kg',
            'U, KJ/KG': 'internal_energy_kj_kg',
            'G, KJ/KG': 'gibbs_kj_kg',
            'S, KJ/(KG)(K)': 'entropy_kj_kg_k',
            'M, (1/n)': 'molecular_weight',
            '(dLV/dLP)t': 'dlv_dlp_t',
            '(dLV/dLT)p': 'dlv_dlt_p',
            'Cp, KJ/(KG)(K)': 'cp_kj_kg_k',
            'GAMMAs': 'gamma',
            'SON VEL,M/SEC': 'sound_speed_m_s',
            'MACH NUMBER': 'mach_number'
        }
        
        for i in range(header_idx + 1, len(lines)):
            line = lines[i].strip()
            if 'PERFORMANCE PARAMETERS' in line or 'MASS FRACTIONS' in line:
                break
                
            # Check if line contains a parameter we care about
            for param_key, param_name in param_mapping.items():
                if line.startswith(param_key):
                    newline = line.replace(param_key, '').strip()
                    values = self._parse_numeric_line(newline)  # Skip parameter name
                    data_dict[param_name] = values 
                    break
        
        # Parse mass fractions
        mass_fractions = self._parse_mass_fractions(section)
        
        # Create states for each location (up to 5 columns)
        num_locations = min(len(locations), len(data_dict.get('pressure_bar', [])))
        
        for j in range(num_locations):
            # Extract mass fractions for this state
            state_mass_fractions = {}
            for species, fractions in mass_fractions.items():
                if j < len(fractions):
                    state_mass_fractions[species] = fractions[j]
            
            state = CEAState(
                location=locations[j],
                pressure_ratio=self._safe_get(data_dict, 'pinf_p_ratio', j, 0.0),
                pressure_bar=self._safe_get(data_dict, 'pressure_bar', j, 0.0),
                temperature_k=self._safe_get(data_dict, 'temperature_k', j, 0.0),
                density_kg_m3=self._parse_scientific(self._safe_get_str(data_dict, 'density_kg_m3', j, '0')),
                enthalpy_kj_kg=self._safe_get(data_dict, 'enthalpy_kj_kg', j, 0.0),
                internal_energy_kj_kg=self._safe_get(data_dict, 'internal_energy_kj_kg', j, 0.0),
                gibbs_kj_kg=self._safe_get(data_dict, 'gibbs_kj_kg', j, 0.0),
                entropy_kj_kg_k=self._safe_get(data_dict, 'entropy_kj_kg_k', j, 0.0),
                molecular_weight=self._safe_get(data_dict, 'molecular_weight', j, 0.0),
                dlv_dlp_t=self._safe_get(data_dict, 'dlv_dlp_t', j, 0.0),
                dlv_dlt_p=self._safe_get(data_dict, 'dlv_dlt_p', j, 0.0),
                cp_kj_kg_k=self._safe_get(data_dict, 'cp_kj_kg_k', j, 0.0),
                gamma=self._safe_get(data_dict, 'gamma', j, 0.0),
                sound_speed_m_s=self._safe_get(data_dict, 'sound_speed_m_s', j, 0.0),
                mach_number=self._safe_get(data_dict, 'mach_number', j, 0.0),
                mass_fractions=state_mass_fractions
            )
            states.append(state)
        
        return states
    
    def _parse_performance(self, section: str) -> List[CEAPerformance]:
        """Parse performance parameters"""
        performance = []
        
        lines = section.split('\n')
        perf_data = {}
        param_map = {
            'Ae/At': 'area_ratio',
            'CSTAR, M/SEC': 'cstar_m_s', 
            'CF': 'cf', 
            'Isp, M/SEC': 'isp_m_s', 
            'Ivac, M/SEC': 'ivac_m_s'
        }
        
        for line in lines:
            line = line.strip()
            for param in param_map.keys():
                if line.startswith(param):
                    
                    values = line.split()
                    if len(values) > 1:
                        #special case: M/SEC can get split due to a space between it and a comma
                        if values[1] == "M/SEC":
                            values[1] = values[0] + " " + values[1]
                            values.pop(0)
                        param_name = param_map[values[0]]
                        perf_data[param_name] = values[1:]
        
        # Create performance objects for each condition
        num_conditions = len(perf_data.get('cstar_m_s', []))
        for i in range(num_conditions):
            perf = CEAPerformance(
                area_ratio=self._safe_get(perf_data, 'area_ratio', i, 0.0),
                cstar_m_s=self._safe_get(perf_data, 'cstar_m_s', i, 0.0),
                cf=self._safe_get(perf_data, 'cf', i, 0.0),
                isp_m_s=self._safe_get(perf_data, 'isp_m_s', i, 0.0),
                ivac_m_s=self._safe_get(perf_data, 'ivac_m_s', i, 0.0)
            )
            performance.append(perf)
        
        return performance
    
    def _parse_mass_fractions(self, section: str) -> Dict[str, List[float]]:
        """Parse mass fractions section"""
        lines = section.split('\n')
        mass_fractions = {}
        
        # Find mass fractions section
        mf_start = -1
        for i, line in enumerate(lines):
            if "MASS FRACTIONS" in line:
                mf_start = i + 2  # Skip header and blank line
                break
        
        if mf_start == -1:
            return mass_fractions
        
        # Parse mass fractions - each line may have multiple species
        for i in range(mf_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('*') and line.endswith('20000.K') or 'NOTE.' in line:
                break
            
            # Split line into parts
            parts = line.split()
            if not parts:
                continue
            
            j = 0
            while j < len(parts):
                # Look for species name (starts with * or is alphabetic)
                if parts[j].startswith('*') or (parts[j].isalpha() and not self._is_numeric(parts[j])) or ('(' in parts[j] and ')' in parts[j]):
                    species = parts[j].replace('*', '')
                    j += 1
                    
                    # Collect numeric values after the species name
                    values = []
                    while j < len(parts) and self._is_numeric(parts[j]):
                        try:
                            values.append(float(parts[j]))
                        except ValueError:
                            break
                        j += 1
                    
                    if values:
                        mass_fractions[species] = values
                else:
                    j += 1
        
        return mass_fractions
    
    def _extract_locations(self, header_line: str) -> List[str]:
        """Extract location names from header"""
        locations = []
        
        # Handle case where header line might not contain location names
        # In CEA output, locations are often implicit from column count
        if "CHAMBER" in header_line:
            locations.append("CHAMBER")
        if "THROAT" in header_line:
            locations.append("THROAT")
        
        # Count EXIT columns
        exit_count = header_line.count("EXIT")
        for i in range(exit_count):
            locations.append(f"EXIT_{i+1}")
        
        # If no explicit location names found, infer from data structure
        # CEA typically has 5 columns: Chamber, Throat, Exit1, Exit2, Exit3
        if not locations:
            locations = ["CHAMBER", "THROAT", "EXIT_1", "EXIT_2", "EXIT_3"]
        
        return locations
    
    def _parse_numeric_line(self, line: str) -> List[str]:
        """Parse a line containing numeric data, preserving scientific notation"""
        # Handle scientific notation like "9.4113 0" -> "9.4113E0"
        result = re.sub(r'(-?\d+\.\d+)\s?([\s-]?\d+)?(?=\s|$)', r'\1E0\2', line)
        result = result.replace('E0-', "E-")
        return result.split()
    
    def _parse_scientific(self, value_str: str) -> float:
        """Parse scientific notation including CEA format like '9.4113 0'"""
        if 'E' in value_str.upper():
            return float(value_str)
        
        # Handle CEA format scientific notation
        parts = value_str.split()
        if len(parts) == 2 and self._is_numeric(parts[0]) and self._is_numeric(parts[1]):
            return float(parts[0]) * (10 ** float(parts[1]))
        
        try:
            return float(value_str)
        except ValueError:
            return 0.0
    
    def _safe_get(self, data_dict: Dict, key: str, index: int, default: float) -> float:
        """Safely get a value from data dictionary"""
        if key not in data_dict or index >= len(data_dict[key]):
            return default
        try:
            value_str = str(data_dict[key][index])
            return self._parse_scientific(value_str)
        except (ValueError, IndexError):
            return default
    
    def _safe_get_str(self, data_dict: Dict, key: str, index: int, default: str) -> str:
        """Safely get a string value from data dictionary"""
        if key not in data_dict or index >= len(data_dict[key]):
            return default
        return str(data_dict[key][index])
    
    def _is_numeric(self, value: str) -> bool:
        """Check if string represents a number"""
        try:
            float(value)
            return True
        except ValueError:
            return False


def main():
    parser = argparse.ArgumentParser(description="Parse CEA output files into JSON test data")
    parser.add_argument("input", type=Path, help="CEA output file to parse")
    parser.add_argument("-o", "--output", type=Path, help="Output JSON file (default: input.json)")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input file {args.input} does not exist")
        return 1
    
    input_path = Path(args.input).resolve()
    output_path = Path(args.output or args.input.with_suffix('.json')).resolve()
    # Parse CEA file
    cea_parser = CEAParser()
    try:
        result = cea_parser.parse_file(args.input)
        
        # Convert to dictionary for JSON serialization
        result_dict = asdict(result)
        
        # Write JSON output
        with open(output_path, 'w') as f:
            if args.pretty:
                json.dump(result_dict, f, indent=2)
            else:
                json.dump(result_dict, f)
        print(f"Successfully parsed {input_path} -> {output_path}")
        # print(f"Found {len(result.equilibrium_states)} equilibrium states")
        # print(f"Found {len(result.frozen_states)} frozen states")
        
    except Exception as e:
        print(f"Error parsing {args.input}: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())