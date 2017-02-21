w = 500
h = 400

svg = d3.select('body').append('svg')
  .attr('width', w)
  .attr('height', h)

#jsonCircles = [
#  { "x_axis": 30, "y_axis": 30, "radius": 20, "color" : "green" },
#  { "x_axis": 70, "y_axis": 70, "radius": 20, "color" : "purple"},
#  { "x_axis": 110, "y_axis": 100, "radius": 20, "color" : "red"}
#]
#
#circles = svg.selectAll("circles")
#  .data(jsonCircles)
#  .enter()
#  .append("circle")
#  .attr("cx", (d) ->  d.x_axis )
#  .attr("cy", (d) ->  d.y_axis )
#  .attr("r", (d) ->  d.radius )
#  .style("fill", (d) ->  d.color )

xScale = d3.scaleLinear()
    .domain([-50, 150])
    .range([50, 400])
#svg.append("g")
#    .call(d3.svg.axis()
#                .scale(xScale)
#                .orient("bottom"))
f = (x, x2, t, y) ->
  svg.append('rect')
      .attr('x', xScale(x))
      .attr('width', xScale(x2) - xScale(x))
      .attr('y', 100)
      .attr('height', 10)
      .attr('opacity', 0.5)
      .attr('fill', 'blue')

  svg.append('text')
      .attr('x', (d) -> xScale((x+x2)/2))
      .attr('y', 100 - y)
      .text(t)
      .attr("text-anchor", "middle")

f(-50, -20, 'Homology', 10)
f(-19, -10, 'long primer', 20)
f(-9, 29, 'extra primer', 30)
f(30, 70, 'two extra primers', 10)
f(71, 125, 'No man', 20)
f(126, 149, 'Synthesis', 30)


