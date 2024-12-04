import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean

class CellTrajectoryAnalyzer:
    def __init__(self):
        self.trajectories = pd.DataFrame()
        self.cell_metadata = pd.DataFrame()
        
    def load_data(self, trajectory_file: str, metadata_file: str = None):
        """
        Load cell trajectory data and optional metadata
        Expected format: frame_number, cell_id, x, y, (optional: z)
        """
        self.trajectories = pd.read_csv(trajectory_file)
        if metadata_file:
            self.cell_metadata = pd.read_csv(metadata_file)
    
    def calculate_motility_metrics(self):
        """Calculate key cell motility metrics"""
        metrics = {}
        
        for cell_id in self.trajectories['cell_id'].unique():
            cell_data = self.trajectories[self.trajectories['cell_id'] == cell_id]
            
            # Sort by frame number
            cell_data = cell_data.sort_values('frame_number')
            
            # Calculate displacement and velocities
            displacements = []
            velocities = []
            
            for i in range(1, len(cell_data)):
                prev = cell_data.iloc[i-1]
                curr = cell_data.iloc[i]
                
                # Calculate displacement
                displacement = euclidean(
                    [prev['x'], prev['y']], 
                    [curr['x'], curr['y']]
                )
                displacements.append(displacement)
                
                # Calculate velocity (assuming constant time between frames)
                velocities.append(displacement)
            
            # Calculate metrics
            metrics[cell_id] = {
                'total_distance': sum(displacements),
                'net_displacement': euclidean(
                    [cell_data.iloc[0]['x'], cell_data.iloc[0]['y']],
                    [cell_data.iloc[-1]['x'], cell_data.iloc[-1]['y']]
                ),
                'average_velocity': np.mean(velocities) if velocities else 0,
                'directional_persistence': self.calculate_persistence(cell_data),
                'track_duration': len(cell_data),
                'confinement_ratio': self.calculate_confinement_ratio(cell_data)
            }
        
        return pd.DataFrame.from_dict(metrics, orient='index')
    
    def calculate_persistence(self, cell_data):
        """Calculate directional persistence"""
        if len(cell_data) < 3:
            return 0
            
        # Calculate angles between consecutive movements
        angles = []
        points = cell_data[['x', 'y']].values
        
        for i in range(len(points)-2):
            v1 = points[i+1] - points[i]
            v2 = points[i+2] - points[i+1]
            
            # Calculate angle between vectors
            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
            angles.append(angle)
            
        return np.mean(np.cos(angles)) if angles else 0
    
    def calculate_confinement_ratio(self, cell_data):
        """Calculate confinement ratio (net displacement / total path length)"""
        points = cell_data[['x', 'y']].values
        
        if len(points) < 2:
            return 0
            
        # Calculate total path length
        total_length = sum(
            np.linalg.norm(points[i+1] - points[i])
            for i in range(len(points)-1)
        )
        
        # Calculate net displacement
        net_displacement = np.linalg.norm(points[-1] - points[0])
        
        return net_displacement / total_length if total_length > 0 else 0
    
    def analyze_population_heterogeneity(self):
        """Analyze heterogeneity in cell population movements"""
        metrics = self.calculate_motility_metrics()
        
        heterogeneity = {
            metric: {
                'mean': metrics[metric].mean(),
                'std': metrics[metric].std(),
                'cv': stats.variation(metrics[metric])  # Coefficient of variation
            }
            for metric in metrics.columns
        }
        
        return pd.DataFrame(heterogeneity)
    
    def plot_trajectories(self, output_file=None):
        """Plot all cell trajectories"""
        plt.figure(figsize=(10, 10))
        
        for cell_id in self.trajectories['cell_id'].unique():
            cell_data = self.trajectories[self.trajectories['cell_id'] == cell_id]
            plt.plot(cell_data['x'], cell_data['y'], '-', alpha=0.5, label=f'Cell {cell_id}')
        
        plt.title('Cell Trajectories')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        
        if output_file:
            plt.savefig(output_file)
        plt.close()
    
    def generate_report(self, output_file="cell_analysis_report.txt"):
        """Generate comprehensive analysis report"""
        metrics = self.calculate_motility_metrics()
        heterogeneity = self.analyze_population_heterogeneity()
        
        with open(output_file, 'w') as f:
            f.write("Cell Trajectory Analysis Report\n")
            f.write("==============================\n\n")
            
            f.write("Population Summary Statistics:\n")
            f.write("-----------------------------\n")
            f.write(heterogeneity.to_string())
            f.write("\n\n")
            
            f.write("Individual Cell Metrics:\n")
            f.write("----------------------\n")
            f.write(metrics.to_string())

# Example usage
if __name__ == "__main__":
    analyzer = CellTrajectoryAnalyzer()
    
    # Load your data
    analyzer.load_data("cell_trajectories.csv", "cell_metadata.csv")
    
    # Generate analysis
    analyzer.generate_report()
    analyzer.plot_trajectories("trajectory_plot.png")
