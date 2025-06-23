# Power System Harmonics Identification Powered by an Eigensystem Realization Approach


## 📄 Manuscript ID:
IEEE LATAM 
## 📬 Submission ID:
XXXX

## ✍️ Authors:
Miguel G. Juarez, Felix Reyes-Maldonado, Alejandro Zamora-Mendez, José Ortiz-Bejar, Juan Carlos Silva-Chavez, Mario R. Arrieta Paternina, and Vicente Torres-Garcia


The novelty of this investigation lies in the power system harmonic identification powered by the ERA method. This approach consists of a system identification technique that allows us to design a filter bank for extracting harmonic components from actual or simulated data, being more flexible and adaptable than the Fourier filters which take into account only a fixed frequency. The ERA filter bank is also able to provide damping and frequency estimates, giving more useful information than Fourier. Thus, the primary contributions of this paper are summarized as follows:


- The estimation achieved by the ERA method allows us to correctly extract harmonics, inter-harmonics, and sub-harmonics from actual or simulated data.

- The proposed method extends the state-of-the-art of parametric/model-based methods by introducing the ERA technique.

- Dynamic modeling: Moves away from static paradigms by employing a dynamic model that captures time-varying harmonic patterns.

- A powerful tool for processing harmonic estimations based on a set of ERA-based FIR filters is proposed to power system community.


## 🛠️ Requirements

- MATLAB: 2018a ✅

## 📂 Project Structure

| Script                      | Related Figure(s)        | Description                                                                                                                   |
|----------------------------|--------------------------|-------------------------------------------------------------------------------------------------------------------------------|
|📁 *Folder*                   | All figures              | *Inrush*, *Theoretical*, *Rectifier*, *WPP_Kundur*: Each folder contains the following files: `ERA_{identifier}.m`, `FFT_{identifier}.m`, `I{identifier}.mat`, and `harmonics{identifier}.mat`. |
| ERA_{identifier}.m         | 1, 4–6, 8–12, 15–24       | Implements the ERA method, reconstructs signals, compares with FFT approach, and plots the Fourier spectrum.                 |
| FFT_{identifier}.m         | 1, 4–6, 8–12, 15–24       | Computes the FFT of the signal and saves the result in `harmonics{identifier}.mat`.                                          |
| I{identifier}.mat          | 2, 4–6, 8, 11, 14, 18, 20, 23 | Stores the original signal for the corresponding study case.                                                              |
| harmonics{identifier}.mat  | 9, 15, 21                 | Stores harmonic identification results related to the study case.                                                            |



## 🧑‍💻 Instructions

### 1.Clone the Repository
Download or clone the Harmonics-LATAM-2025 repository:

``` git clone https://github.com/your-username/Harmonics-LATAM-2025.git```

### 2.Explore Study Cases
Inside the repository, you'll find four folders, each corresponding to a different case study presented in the article.
Choose the one that best fits your area of interest.

### 3.Run the ERA Algorithm
In the selected folder, locate the file named ERA.m and execute it using MATLAB.
This script applies the ERA-based filter bank to extract and analyze harmonic components.

### 3.Analyze Results
Review the simulation outputs and plots corresponding to the chosen study case to see how the method performs on real or simulated data.
