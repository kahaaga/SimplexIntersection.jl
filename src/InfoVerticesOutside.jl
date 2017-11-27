function InfoVerticesOutside(SimplexVert,SimplexVertTarget,IndexNotContainedVertices,SimplexOrientation,tolerance)
  #SimplexVert : Matrix of vertices of the holding simplex
  #nx(n+1) matrix containing the set of vertices furnishing the simplex. The vertices are the columns

  #SimplexVertTarget : Matrix of vertices of the simplex which conves expansion in terms of the vertices of SimplexVert we want to obtain

  #SimplexOrientation = determinat of the vectors generating the simplex in the
  #order given by SimplexVert and taking the first one as the origin

  #IndexNotContainedVertices : Row vector of indices of vertices in SimplexVertTarget lying outside the circumsphere of SimplexVert

  #Outcome:

  # InfoVertices is a matrix of dimension (n+2) x size(SimplexVertTarget,2)

  #1st row corresponds to the indices of the analized vertices.

  #following rows:

  #Below each index there is a column with the indices of the beta coefficients in its convex expansion

  #Example in 3d: say that the vertex of S1 with index 1 belongs has the
  #expansion (v1)_1 = beta(1)_1 (v2)_1 + beta(1)_2 (v2)_2 + beta(1)_3 (v2)_3 + beta(1)_4 (v2)_4
  #in the vertices of S2

  #The the result in this case would read

  #[1;beta(1)_1;beta(1)_2;beta(1)_3;beta(1)_4]


  #Notice that, in this case, the sign of beta(i)_j says to which side of the j-th face of S2, is the i-th vertex of S1 lying, according to the orientaion inherited by S2

  n = size(SimplexVert, 1)
  Points = copy(SimplexVertTarget[:, IndexNotContainedVertices])
  np = size(Points, 2)
  InfoVertices = zeros(n + 2, np)
  InfoVertices[1, :] = IndexNotContainedVertices
  for i = 1:np
    for j = 1:(n + 1)
      Aux = copy(SimplexVert)
      Aux[:, j] = Points[:, i]
      beta = det([ones(1, n + 1); Aux])/SimplexOrientation
      if abs(beta) <= tolerance
        beta = 0
      end
      InfoVertices[j + 1, i] = beta
    end
  end
  return InfoVertices
end
