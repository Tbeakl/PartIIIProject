# using MultivariateStats, RDatasets

# iris = dataset("datasets", "iris")

# X = Matrix(iris[1:2:end,1:4])'
# X_labels = Vector(iris[1:2:end,5])

# pca = fit(PCA, X)
# println("test")
# Ypca = predict(pca, X)
using DSP, FFTW

c_in = [0.5, 0.5, 0., 0., 0., 0., 0., 0.]
a = [0.25, 0.25, 0.25, 0.25, 0., 0., 0., 0.]
b = [0.25, 0.25, 0.25, 0.25, 0., 0., 0., 0.]
out = [.125, .125, .125, .125, .125, .125, .125, .125]

a_tild = fft(a)
b_tild = fft(b)
c_tild = fft(c_in)
out_tild = fft(out)

t_a = real(ifft(conj.(b_tild) .* conj.(c_tild) .* out_tild))
t_b = real(ifft(conj.(a_tild) .* conj.(c_tild) .* out_tild))
t_c = real(ifft(conj.(b_tild .* a_tild) .* out_tild))
t_out = real(ifft(b_tild .* c_tild .* a_tild))