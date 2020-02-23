from patch_util.patch import design_patch, input_impedance, inset_feed_position, \
    get_directivity, patch_eh_plane_plot, surface_plot, getGs


freq = 2.4e9
Er = 4.4
h = 1.6 * 10 ** -3
v = 3 * 10 ** 8

W, L = design_patch(Er, h, freq)


Rin = input_impedance(freq, W, L)
print('Inset Feed Position : ', inset_feed_position(Rin, L))

G1, G12 = getGs(freq, W, L)
print('G1 : ', G1)
print('G12 : ', G12)

I1 = 1.863
I2 = 3.59801

d1, d2 = get_directivity(G1, G12, W, freq, I1, I2)

print('Directivity : ', d1, ' dB')
print('Directivity : ', d2, ' dB')

fields = patch_eh_plane_plot(freq, W, L, h, Er)
surface_plot(fields)
