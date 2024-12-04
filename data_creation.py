import pandas as pd
import numpy as np
from datetime import datetime, timedelta

def generate_trajectory_data(num_people=100, points_per_trajectory=50, bias_factors=True):
    """
    Generate synthetic trajectory data with optional built-in bias factors.
    
    Parameters:
    num_people: Number of unique individuals
    points_per_trajectory: Number of points per person's trajectory
    bias_factors: If True, introduce systematic biases for testing
    """
    # Create lists to store data
    data = []
    
    for person_id in range(num_people):
        # Generate true trajectory
        t = np.linspace(0, 10, points_per_trajectory)
        
        # True path is a 3D curve with some randomness
        x = 2 * np.cos(t) + np.random.normal(0, 0.1, points_per_trajectory)
        y = 2 * np.sin(t) + np.random.normal(0, 0.1, points_per_trajectory)
        z = t/3 + np.random.normal(0, 0.1, points_per_trajectory)
        
        # Generate predictions with intentional biases if enabled
        if bias_factors and person_id < num_people // 2:
            # Introduce systematic bias for first half of people
            predicted_x = x + 0.2 + np.random.normal(0, 0.2, points_per_trajectory)
            predicted_y = y - 0.1 + np.random.normal(0, 0.2, points_per_trajectory)
            predicted_z = z + 0.3 + np.random.normal(0, 0.2, points_per_trajectory)
        else:
            # Less bias for second half
            predicted_x = x + np.random.normal(0, 0.1, points_per_trajectory)
            predicted_y = y + np.random.normal(0, 0.1, points_per_trajectory)
            predicted_z = z + np.random.normal(0, 0.1, points_per_trajectory)
        
        # Generate timestamps
        start_time = datetime(2024, 1, 1)
        timestamps = [start_time + timedelta(seconds=i*30) for i in range(points_per_trajectory)]
        
        # Combine into data points
        for i in range(points_per_trajectory):
            data.append({
                'person_id': person_id,
                'timestamp': timestamps[i],
                'x': x[i],
                'y': y[i],
                'z': z[i],
                'predicted_x': predicted_x[i],
                'predicted_y': predicted_y[i],
                'predicted_z': predicted_z[i]
            })
    
    return pd.DataFrame(data)

def generate_demographic_data(num_people=100, bias_factors=True):
    """
    Generate synthetic demographic data with optional built-in bias factors.
    """
    # Define possible values for each demographic category
    age_groups = ['18-25', '26-35', '36-45', '46-55', '56+']
    genders = ['Male', 'Female', 'Non-binary']
    ethnicities = ['Group A', 'Group B', 'Group C', 'Group D']
    
    # Generate random demographic data
    data = []
    for person_id in range(num_people):
        # Introduce some correlations if bias_factors is True
        if bias_factors and person_id < num_people // 2:
            # First half of people more likely to be in certain groups
            age_group = np.random.choice(age_groups[:2])  # Younger age groups
            gender = np.random.choice(genders[:2])  # Binary genders
            ethnicity = np.random.choice(ethnicities[:2])  # First two ethnic groups
        else:
            # Second half more evenly distributed
            age_group = np.random.choice(age_groups)
            gender = np.random.choice(genders)
            ethnicity = np.random.choice(ethnicities)
        
        data.append({
            'person_id': person_id,
            'age_group': age_group,
            'gender': gender,
            'ethnicity': ethnicity
        })
    
    return pd.DataFrame(data)

# Generate and save the data
if __name__ == "__main__":
    num_people = 100
    
    # Generate trajectory data
    trajectory_df = generate_trajectory_data(num_people=num_people, bias_factors=True)
    trajectory_df.to_csv('trajectory_data.csv', index=False)
    
    # Generate demographic data
    demographic_df = generate_demographic_data(num_people=num_people, bias_factors=True)
    demographic_df.to_csv('demographic_data.csv', index=False)
    
    print("Generated data files:")
    print("- trajectory_data.csv")
    print("- demographic_data.csv")
