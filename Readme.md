# Cantilever Beam Topology Optimization

This MATLAB code performs topology optimization for a cantilever beam using the level set method. It aims to minimize compliance while maintaining a target volume fraction.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Boundary Conditions](#boundary-conditions)
- [Running the Optimization](#running-the-optimization)
- [Visualization](#visualization)
- [Output](#output)
- [Customization](#customization)
- [License](#license)

## Installation

1. Ensure you have MATLAB 2022b or higher installed on your system.
2. Clone this repository.

```bash
git clone https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness.git
```

## Usage

1. Open MATLAB and navigate to the directory containing `Canteliver_Moment.m`.
2. Run the script:

```matlab
run Canteliver_Moment.m
```

## Configuration

Adjust the main parameters at the beginning of the script:

```matlab
MP = struct("Emax",1,"Emin",1e-3,"v",0.3);
GP = struct("ngp",4,"nor",1,"dx",1,"dy",1,"X0",0,"Y0",30,"lx",60,"ly",30);
OPP = struct("tho",0.5,"Noi",1000,"dt",0.05,"SP",1);
```

- `MP`: Material properties
- `GP`: Geometry parameters
- `OPP`: Optimization parameters

## Boundary Conditions

Modify the Essential Boundary Conditions (EBC) and Natural Boundary Conditions (NBC) in the script:

```matlab
indEBC = find(NodeCoord(:,1)==0);
EBC = [indEBC,zeros(numel(indEBC),2)];
indNBC = find(abs(NodeCoord(:,2)-15)<=0.1 & NodeCoord(:,1)==60);
NBC = [indNBC,zeros(numel(indNBC),1),-1/numel(indNBC)*ones(numel(indNBC),1)];
```

## Running the Optimization

1. Set the desired number of iterations in `OPP.Noi`.
2. Run the script in MATLAB.
3. The optimization will proceed, displaying progress at each iteration.

## Visualization

Set `OPP.SP = 1` to display the evolving optimized shape after each iteration.
|Initial design | Optimum shape |
| ------------- | ------------- |
|Full Solid     |![Fullsolid](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/498552ef-2236-4742-b057-e5ac1f12803a)|
|Case1          | ![Case1](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/b299eb36-4c52-476e-9657-7722aa4e01cf)|
|Case2          | ![Case2](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/a1c4837e-4fa6-4f58-ad66-ddb91505a1cf)|
|Case3          | ![Case3](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/f0b8bbac-a963-49e5-9840-be6d90d40278)|
|Case4          | ![Case4](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/4d079a8a-9d22-4740-9b63-7195cdd8b216)|
|Case5          | ![Case5](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/3a1ffe9b-44e5-4f8d-b6ca-2a3cca2b8026)|
|Case6          | ![Case6](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/51e3e730-ca88-49ce-b2a1-59464fbd1c61)|
|Case7          | ![Case7](https://github.com/msdamzdh/TOLSF_InitialDesignEffectiveness/assets/155764233/8166d85e-a52b-46ab-bffc-03cf716f847a)|

## Output

The script outputs:
- Optimized shape of the cantilever beam
- Convergence history of the objective function and volume fraction

## Customization

You can modify various aspects of the code:
- Mesh resolution: Change `dx` and `dy`
- Domain size: Adjust `lx` and `ly`
- Material properties: Modify the `MP` structure
- Optimization parameters: Update the `OPP` structure

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
