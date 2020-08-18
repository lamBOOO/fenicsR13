import numpy

dofs = 1425
kn_range1 = [1.0]
numknudsen1 = len(kn_range1)
kn_range2 = [(i + 1) / 10 for i in range(10)]
numknudsen2 = len(kn_range2)
fields = [
    "theta", "s_x", "s_y", "p", "u_x", "u_x", "sigma_xx", "sigma_xy", "sigma_yx"
]
numfields = len(fields)

V = numpy.zeros(shape=(numknudsen1 * numknudsen2, numfields * dofs))

for iKn, kn1 in enumerate(kn_range1):
    for jKn, kn2 in enumerate(kn_range2):
        # row = numpy.zeros(shape=(numfields * dofs))

        field_data = numpy.zeros(shape=(numfields, dofs))
        for i, field in enumerate(fields):
            data = numpy.loadtxt("{}_{}/{}_0.mat".format(kn1, kn2, field))
            field_data[i] = data
            # print(field_data)
            # print(data)
            # print(len(data))
            # print(data.shape)
            # print("next")
        row = field_data.reshape(numfields * dofs)
        # print(row)
        print(iKn * numknudsen1 + jKn)
        V[iKn * numknudsen1 + jKn] = row

numpy.savetxt("V_{}.mat".format(numknudsen1 * numknudsen2), V)
print(V)
