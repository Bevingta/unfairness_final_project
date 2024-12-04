import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy import stats
import seaborn as sns
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt

class TrajectoryFairnessAnalyzer:
    def __init__(self):
        self.trajectories = pd.DataFrame()
        self.demographic_data = pd.DataFrame()
        self.error_metrics = {}
        
    def load_data(self, trajectory_file: str, demographic_file: str) -> None:
        """
        Load trajectory and demographic data from files.
        
        Parameters:
        trajectory_file: Path to CSV with columns [person_id, timestamp, x, y, z, predicted_x, predicted_y, predicted_z]
        demographic_file: Path to CSV with columns [person_id, age_group, gender, ethnicity, etc.]
        """
        self.trajectories = pd.read_csv(trajectory_file)
        self.demographic_data = pd.read_csv(demographic_file)
        
    def calculate_errors(self) -> Dict:
        """Calculate trajectory prediction errors for each person"""
        errors = {}
        
        for person_id in self.trajectories['person_id'].unique():
            person_data = self.trajectories[self.trajectories['person_id'] == person_id]
            
            # Calculate various error metrics
            errors[person_id] = {
                'mse_x': mean_squared_error(person_data['x'], person_data['predicted_x']),
                'mse_y': mean_squared_error(person_data['y'], person_data['predicted_y']),
                'mse_z': mean_squared_error(person_data['z'], person_data['predicted_z']),
                'mae_x': mean_absolute_error(person_data['x'], person_data['predicted_x']),
                'mae_y': mean_absolute_error(person_data['y'], person_data['predicted_y']),
                'mae_z': mean_absolute_error(person_data['z'], person_data['predicted_z']),
                'max_error': np.max(np.abs(np.array([
                    person_data['x'] - person_data['predicted_x'],
                    person_data['y'] - person_data['predicted_y'],
                    person_data['z'] - person_data['predicted_z']
                ])))
            }
            
        return errors

    def analyze_bias(self) -> pd.DataFrame:
        """Analyze error distributions across demographic groups"""
        errors = self.calculate_errors()
        
        # Create DataFrame with errors and demographic info
        error_df = pd.DataFrame.from_dict(errors, orient='index')
        analysis_df = error_df.merge(self.demographic_data, left_index=True, right_on='person_id')
        
        # Initialize results dictionary
        bias_analysis = {}
        
        # Analyze each demographic category
        demographic_cols = ['age_group', 'gender', 'ethnicity']  # Add more as needed
        error_metrics = ['mse_x', 'mse_y', 'mse_z', 'mae_x', 'mae_y', 'mae_z', 'max_error']
        
        for demo_col in demographic_cols:
            for error_metric in error_metrics:
                # Calculate summary statistics for each group
                group_stats = analysis_df.groupby(demo_col)[error_metric].agg([
                    'mean', 'std', 'median', 'count'
                ])
                
                # Perform statistical tests
                groups = [group for _, group in analysis_df.groupby(demo_col)[error_metric]]
                f_stat, p_value = stats.f_oneway(*groups)
                
                bias_analysis[f"{demo_col}_{error_metric}"] = {
                    'stats': group_stats,
                    'f_statistic': f_stat,
                    'p_value': p_value
                }
        
        return pd.DataFrame(bias_analysis)

    def identify_common_mistakes(self) -> Dict:
        """Identify patterns in prediction errors across different groups"""
        mistakes = {}
        
        for demo_col in ['age_group', 'gender', 'ethnicity']:
            group_mistakes = {}
            
            for group in self.demographic_data[demo_col].unique():
                # Get person_ids for this group
                group_ids = self.demographic_data[
                    self.demographic_data[demo_col] == group
                ]['person_id']
                
                # Get trajectory data for this group
                group_trajectories = self.trajectories[
                    self.trajectories['person_id'].isin(group_ids)
                ]
                
                # Calculate error patterns
                error_patterns = {
                    'systematic_bias_x': np.mean(
                        group_trajectories['predicted_x'] - group_trajectories['x']
                    ),
                    'systematic_bias_y': np.mean(
                        group_trajectories['predicted_y'] - group_trajectories['y']
                    ),
                    'systematic_bias_z': np.mean(
                        group_trajectories['predicted_z'] - group_trajectories['z']
                    ),
                    'variance_x': np.var(
                        group_trajectories['predicted_x'] - group_trajectories['x']
                    ),
                    'variance_y': np.var(
                        group_trajectories['predicted_y'] - group_trajectories['y']
                    ),
                    'variance_z': np.var(
                        group_trajectories['predicted_z'] - group_trajectories['z']
                    )
                }
                
                group_mistakes[group] = error_patterns
            
            mistakes[demo_col] = group_mistakes
            
        return mistakes

    def generate_report(self, output_file: str = "fairness_report.txt") -> None:
        """Generate a comprehensive fairness analysis report"""
        bias_analysis = self.analyze_bias()
        common_mistakes = self.identify_common_mistakes()
        
        with open(output_file, 'w') as f:
            f.write("Trajectory Analysis Fairness Report\n")
            f.write("==================================\n\n")
            
            # Write bias analysis results
            f.write("Bias Analysis Results:\n")
            f.write("---------------------\n")
            for col in bias_analysis.columns:
                f.write(f"\n{col}:\n")
                f.write(bias_analysis[col].to_string())
                f.write("\n")
            
            # Write common mistakes analysis
            f.write("\nCommon Mistakes Analysis:\n")
            f.write("------------------------\n")
            for demo_col, mistakes in common_mistakes.items():
                f.write(f"\n{demo_col}:\n")
                for group, patterns in mistakes.items():
                    f.write(f"\n  {group}:\n")
                    for pattern, value in patterns.items():
                        f.write(f"    {pattern}: {value:.4f}\n")

# Example usage
if __name__ == "__main__":
    analyzer = TrajectoryFairnessAnalyzer()
    
    # Load your data
    analyzer.load_data(
        trajectory_file="trajectory_data.csv",
        demographic_file="demographic_data.csv"
    )
    
    # Generate comprehensive analysis
    analyzer.generate_report()
