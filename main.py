# By Lishen Wang, May 16 2023, reference to https://www.mathworks.com/matlabcentral/fileexchange/106870-calculate-langlands-automorphic-mass-charge-spectrum

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def Langlands(Wuhan, z):
    # Automorphic function
    Seq = (Wuhan[:, 2] + z * Wuhan[:, 3]) / (Wuhan[:, 0] + z * Wuhan[:, 1])
    # Average of mass enthalpy to charge potential
    avg = np.mean(Seq)
    # Fractal moment for 11th string world
    outputArg1 = np.sum((Seq - avg) ** (1 / 11)) / (len(Seq) - 1)
    return complex(outputArg1)


def getVirusFourLabel(oneVirusSeq):
    # Charge, pH, Mass, Hydro indice
    oneVirusSeqE = oneVirusSeq[:, 0]
    oneVirusSeqPH = oneVirusSeq[:, 1]
    oneVirusSeqM = oneVirusSeq[:, 2]
    oneVirusSeqHydro = oneVirusSeq[:, 3]

    # Remove NaN values
    oneVirusSeqE = oneVirusSeqE[~np.isnan(oneVirusSeqE)]
    oneVirusSeqPH = oneVirusSeqPH[~np.isnan(oneVirusSeqPH)]
    oneVirusSeqM = oneVirusSeqM[~np.isnan(oneVirusSeqM)]
    oneVirusSeqHydro = oneVirusSeqHydro[~np.isnan(oneVirusSeqHydro)]

    oneVirus = np.column_stack((oneVirusSeqE, oneVirusSeqPH, oneVirusSeqM, oneVirusSeqHydro))

    return oneVirus


# read data from Excelsheet
data = pd.read_excel('NewSpikeAnalysisV17murineOC43W2Qv2.xls', sheet_name=2, header=1).values

nameVirus = ['Omicron', 'Wuhan', 'Delta', 'XXXX', 'France', 'XXXXX', 'Murine', 'XXXXX', 'OC43']
numAcid, numVirus = data.shape
numVirus = numVirus // 5

for i in range(int(9)):
    # extract the data
    colStart = i * 5
    colEnd = i * 5 + 4
    oneVirus = data[:, colStart:colEnd]
    nameV = nameVirus[i]
    exec(f"{nameV} = oneVirus")

    del oneVirus

value = 0
# Divide into s1 and s2 domains roughly
OmicronFront = Omicron[0:len(Omicron) // 2 - value, :]
WuhanFront = Wuhan[0:len(Wuhan) // 2 - value, :]
DeltaFront = Delta[0:len(Delta) // 2 - value, :]
OC43Front = OC43[0:len(OC43) // 2 - value, :]
MurineFront = Murine[0:len(Murine) // 2 - value, :]

OmicronEnd = Omicron[len(Omicron) // 2 - value:, :]
WuhanEnd = Wuhan[len(Wuhan) // 2 - value:, :]
DeltaEnd = Delta[len(Delta) // 2 - value:, :]
OC43End = OC43[len(OC43) // 2 - value:, :]
MurineEnd = Murine[len(Murine) // 2 - value:, :]

# Run equivalent receptor-binding domain distance from 0.1 to 10 nm
# with 1pm step for s1 domain
recordZ = np.arange(0.1, 10.001, 0.001)
record = np.zeros(len(recordZ), dtype=complex)
record2 = np.zeros(len(recordZ), dtype=complex)
record3 = np.zeros(len(recordZ), dtype=complex)
record4 = np.zeros(len(recordZ), dtype=complex)
record5 = np.zeros(len(recordZ), dtype=complex)

for i, z in enumerate(recordZ):
    record[i] = np.log(Langlands(OmicronFront, z))
    record2[i] = np.log(Langlands(WuhanFront, z))
    record3[i] = np.log(Langlands(DeltaFront, z))
    record4[i] = np.log(Langlands(OC43Front, z))
    record5[i] = np.log(Langlands(MurineFront, z))
print(record[1])

# The magnitude of spectrum
plt.figure(1)
plt.loglog(recordZ, np.abs(record), recordZ, np.abs(record2), '--', recordZ, np.abs(record3), recordZ, np.abs(record4),
           '--', recordZ, np.abs(record5), linewidth=1)
plt.legend(['Omicron', 'Wuhan', 'Delta', 'OC43', 'Murine'])
plt.gcf().set_size_inches(13, 3)
plt.gca().set_position([0.04, 0.14, 0.94, 0.8])
plt.gca().tick_params(axis='x')
plt.gca().tick_params(axis='y')
plt.gca().tick_params(labelsize=16)
plt.savefig('Fig 1 mag effects of different Z of the Langlands program on viruses', dpi=600)
plt.show()

# The phase of spectrum
plt.figure(2)
plt.loglog(recordZ, np.angle(record), recordZ, np.angle(record2), '--', recordZ, np.angle(record3), recordZ,
           np.angle(record4), '--', recordZ, np.angle(record5), linewidth=1)
plt.legend(['Omicron', 'Wuhan', 'Delta', 'OC43', 'Murine'])
plt.gcf().set_size_inches(13, 3)
plt.gca().set_position([0.04, 0.14, 0.94, 0.8])
plt.gca().tick_params(axis='x')
plt.gca().tick_params(axis='y')
plt.gca().tick_params(labelsize=16)
plt.savefig('Fig 2 ang effects of different Z of the Langlands program on viruses', dpi=600)
plt.show()
