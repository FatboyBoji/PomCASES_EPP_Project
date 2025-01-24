import numpy as np

# Base price ranges (€/kg)
green_h2_price = np.random.uniform(4, 8, 24)  # Green hydrogen
blue_h2_price = np.random.uniform(2, 4, 24)   # Blue hydrogen

# Introduce fluctuations
fluctuations = np.sin(np.linspace(0, 2 * np.pi, 24)) * 0.5  # Daily cycle

# Construct test data
c_hBuy = green_h2_price + fluctuations
c_hSell = c_hBuy * 0.6  # Assume selling price is 60% of buying price

print("c_hBuy:", c_hBuy)
print("c_hSell:", c_hSell)
