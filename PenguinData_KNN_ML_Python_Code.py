# Cody Helm

# Import libraries
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.metrics import accuracy_score

# Import dataset
df = pd.read_csv('https://raw.githubusercontent.com/mwaskom/seaborn-data/master/penguins.csv')
df.head()


#%%
# Delete sex variable
#del df['sex']

# Drop missing rows
df.dropna(axis = 0, how = 'any', subset = None, inplace = True)

# Convert non-numeric data using one-hot encoding
df = pd.get_dummies(df, columns=['island','sex'])

# Scale the independent variables and ignore the variable/column (sex)
scaler = StandardScaler()
scaler.fit(df.drop('species',axis=1))
scaled_df = scaler.transform(df.drop('species',axis=1))

# Assign X and y variables
X = scaled_df
y = df['species'] # species is what we are trying to guess here


# %%
# Split data into test/train set (70/30 split) and shuffle
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, shuffle=True)

# Assign KNeighborsClassifier algorithm
my_model = KNeighborsClassifier(n_neighbors=5)

# Fit algorithm to training dataset
my_model.fit(X_train, y_train)

# Run algorithm on test dataset
my_model_test = my_model.predict(X_test)


#%% 
# Evaluate predictions

accuracy = accuracy_score(y_test, my_model_test)
print(f"Model Accuracy: {accuracy:.4f}")

print(f"Confusion Matrix: ", confusion_matrix(y_test,my_model_test))
print(f"Classification Report: \n", classification_report(y_test,my_model_test))


#%%
# Data point to predict
penguin = pd.DataFrame([[
    39, #bill_length
    18, #bill_depth
    190, #flipper_length
    3500, #body_mass_g
    0, #island_Biscoe
    0, #island_Dream
    1, #island_Torgersen
    0, #sex_Female
    1, #sex_Male
]],columns=['bill_length_mm', 'bill_depth_mm', 
            'flipper_length_mm', 'body_mass_g','island_Biscoe','island_Dream',
            'island_Torgersen','sex_FEMALE','sex_MALE']
)

# penguin = pd.DataFrame([[39, 18, 190, 3500, 0, 0, 1, 0, 1]])

penguin_scaled = scaler.transform(penguin)

# Make prediction
predict_penguin = my_model.predict(penguin_scaled)
predict_penguin

print(f"Prediction for new penguin is: {predict_penguin} and the true penguin is Adelie.")