# Matlab_CVA_DVA_Calculation
This MATLAB program for performing Credit Valuation Adjustment (CVA) and Debit Valuation Adjustment (DVA) calculations, including market simulation, portfolio pricing, and exposure assessment.


## Installation
1. Ensure MATLAB is installed (tested on MATLAB R2020a or newer).
2. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/CVA_Calculation_MATLAB.git
   ```
3. Add the repository folder and subfolders to MATLAB's search path:
   ```matlab
   addpath(genpath('CVA_Calculation_MATLAB'));
   ```

## Usage
### Running the Program
1. Launch MATLAB.
2. Execute the `Main.m` script.
3. When prompted, select the input data file containing:
   - Market data (e.g., yield curves, credit spreads).
   - Portfolio details (e.g., forward and swap transactions).
4. The program will:
   - Simulate market scenarios.
   - Perform portfolio pricing and exposure calculations.
   - Output results to an Excel file.

### Input and Output
- **Input**: Raw data file containing market and portfolio information.
- **Output**: Excel file summarizing:
  - Netting set-level metrics: EPE, PFE, CVA, ENE, PNE, and DVA.
  - Deal-level metrics: MTM, EPE, PFE, CVA, ENE, PNE, and DVA.

## Business Context
The program assists financial institutions in:
- Evaluating counterparty credit risk (CVA/DVA).
- Analyzing netting benefits from legal agreements.
- Reporting key exposure metrics for regulatory compliance (e.g., Basel requirements).

### Detailed Sections
#### Market Environment Development
- Sets up the valuation date and market conditions.
- Ensures consistent calculations for pricing and risk analysis.

#### Portfolio Data
- Reads transactions and counterparty data.
- Establishes the portfolio structure for CVA/DVA analysis.

#### Simulation
- Generates market scenarios for future price and exposure predictions.
- Supports risk assessment under varying conditions.

#### Exposure and CVA/DVA Calculation
- Values financial instruments based on simulated market data.
- Computes CVA/DVA and aggregates metrics at the netting set level.
- Allocates risk metrics to individual deals.


