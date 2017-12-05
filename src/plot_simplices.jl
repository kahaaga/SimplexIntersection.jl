using PlotlyJS

#' Plots a pair of simplices given as dim x n_vertices Float64 arrays.
function plot_simplex(simplex)

   data = PlotlyJS.GenericTrace[]

   layout = PlotlyJS.Layout(autosize = true,
                             margin = PlotlyJS.attr(l=0, r=0, b=0, t=65))


   s1 = hcat(simplex, simplex[1:end, 2])
   s1 = hcat(s1, s1[1:end, 3])
   s1 = hcat(s1, s1[1:end, 1])
   s1 = hcat(s1, s1[1:end, 4])

   trace = PlotlyJS.scatter3d(;x = s1[1, 1:end],
                              y = s1[2, 1:end],
                              z = s1[3, 1:end],
                              mode = "markers+lines",
                              showlegend = false,
                              marker = PlotlyJS.attr(size = 3, color = "black"),
                              line = PlotlyJS.attr(size = 1, color = "white"))

   push!(data, trace)

   PlotlyJS.plot(data, layout)
end



#' Plots a pair of simplices given as dim x n_vertices Float64 arrays.
function plot_simplices(simplex1, simplex2)

   data = PlotlyJS.GenericTrace[]

   layout = PlotlyJS.Layout(autosize = true,
                             margin = PlotlyJS.attr(l=0, r=0, b=0, t=65))


   s1 = hcat(simplex1, simplex1[1:end, 2])
   s1 = hcat(s1, s1[1:end, 3])
   s1 = hcat(s1, s1[1:end, 1])
   s1 = hcat(s1, s1[1:end, 4])

   trace = PlotlyJS.scatter3d(;x = s1[1, 1:end],
                              y = s1[2, 1:end],
                              z = s1[3, 1:end],
                              mode = "markers+lines",
                              showlegend = false,
                              marker = PlotlyJS.attr(size = 3, color = "black"),
                              line = PlotlyJS.attr(size = 1, color = "white"))

   push!(data, trace)

   s2 = hcat(simplex2, simplex2[1:end, 2])
   s2 = hcat(s2, s2[1:end, 3])
   s2 = hcat(s2, s2[1:end, 1])
   s2 = hcat(s2, s2[1:end, 4])

   trace = PlotlyJS.scatter3d(;x = s2[1, 1:end],
                              y = s2[2, 1:end],
                              z = s2[3, 1:end],
                              mode = "markers+lines",
                              showlegend = false,
                              marker = PlotlyJS.attr(size = 3, color = "black"),
                              line = PlotlyJS.attr(size = 1, color = "blue"))

   push!(data, trace)

   PlotlyJS.plot(data, layout)
end



#' Plots a pair of simplices given as dim x n_vertices Float64 arrays.
function plot_simplex_and_triangulation(simplex1, triang_vertices, triang_simplex_indices)

   data = PlotlyJS.GenericTrace[]

   layout = PlotlyJS.Layout(autosize = true,
                             margin = PlotlyJS.attr(l=0, r=0, b=0, t=65))


   s1 = hcat(simplex1, simplex1[1:end, 2])
   s1 = hcat(s1, s1[1:end, 3])
   s1 = hcat(s1, s1[1:end, 1])
   s1 = hcat(s1, s1[1:end, 4])

   trace = PlotlyJS.scatter3d(;x = s1[1, 1:end],
                              y = s1[2, 1:end],
                              z = s1[3, 1:end],
                              mode = "markers+lines",
                              showlegend = false,
                              marker = PlotlyJS.attr(size = 3, color = "blue"),
                              line = PlotlyJS.attr(size = 1, color = "blue"))

   push!(data, trace)

   for i = 1:size(triang_simplex_indices, 1)
      s = triang_vertices[triang_simplex_indices[i, 1:end], 1:end].'
      s = hcat(s, s[1:end, 2])
      s = hcat(s, s[1:end, 3])
      s = hcat(s, s[1:end, 1])
      s = hcat(s, s[1:end, 4])


      trace = PlotlyJS.scatter3d(;x = s[1, 1:end],
                           y = s[2, 1:end],
                           z = s[3, 1:end],
                           mode = "markers+lines",
                           showlegend = false,
                           marker = PlotlyJS.attr(size = 3, color = "black"),
                           line = PlotlyJS.attr(size = 1, color = "white"))

      push!(data, trace)
   end


   PlotlyJS.plot(data, layout)
end
