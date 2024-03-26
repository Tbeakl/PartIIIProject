using MultivariateStats, RDatasets

iris = dataset("datasets", "iris")

X = Matrix(iris[1:2:end,1:4])'
X_labels = Vector(iris[1:2:end,5])

pca = fit(PCA, X)
println("test")
Ypca = predict(pca, X)