w = 960
h = 5000

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

xpadding = 50.0


tooltip = d3.select("body").append("div")
  .attr("class", "tooltip")


d3.json('data.json', (data) ->

  query_len = data.contigs[0].query_length

  xScale = d3.scaleLinear()
    .domain([0, query_len])
    .range([0+xpadding, w-xpadding])

  yScale = d3.scaleLinear()
    .domain([0, data.contigs.length])
    .range([120, 500])

  svg.append('rect')
    .attr('x', xScale(0))
    .attr('width', xScale(query_len))
    .attr('y', 100)
    .attr('height', 10)
    .attr('opacity', 0.5)
    .attr('fill', 'blue')

  svg.selectAll('rects')
    .data(data.contigs)
    .enter()
    .append('rect')
    .attr('x', (d) -> xScale(d.q_start))
    .attr('y', (d,i) -> 120 + i * 10)
    .attr('width', (d) -> xScale(d.q_end) - xScale(d.q_start))
    .attr('height', 7)
    .attr('fill', fill)
    .attr('opacity', 0.9)
    .on("mouseover", (d) ->
      coordinates = d3.mouse(this)
      d3.select(this).style("fill", 'red')
      tooltip.text(d.subject_acc + '<br>' + d.q_start + ', ' + d.q_end + '<br>' + d.contig_type)
        .style("visibility", "visible")
      )
    .on("mousemove", () ->
      tooltip.style("top", (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px"))
    .on("mouseout", (d) ->
      d3.select(this).style('fill', fill)
      tooltip.style("visibility", "hidden")
  )

  svg.selectAll('rects')
    .data(data.contigs)
    .enter()
    .append('rect')
    .attr('x', (d) -> xScale(d.q_start))
    .attr('y', (d,i) -> 100)
    .attr('width', 1)
    .attr('height', h)
    .attr('fill', 'black')
    .attr('opacity', 0.1)

  svg.selectAll('rects')
    .data(data.contigs)
    .enter()
    .append('rect')
    .attr('x', (d) -> xScale(d.q_end))
    .attr('y', (d,i) -> 100)
    .attr('width', 1)
    .attr('height', h)
    .attr('fill', 'blue')
    .attr('opacity', 0.1)

  d3.json('primer_data.json', (primerdata) ->
      svg.selectAll('primer_rects')
        .data(primerdata.contigs)
        .enter()
        .append('rect')
        .attr('x', (d) -> xScale(d.q_start))
        .attr('y', (d,i) -> 100)
        .attr('width', 1)
        .attr('height', h)
        .attr('fill', "red")
        .attr('opacity', 0.35)
  )

)



fill = (d) ->
  return 'black' if d.contig_type == 'contig'
  return 'purple' if d.contig_type == 'product'
  return 'blue' if d.contig_type == 'gap'
  return 'orange'


#width = 200
#tspan = text.text(null).append('tspan')
#tspan.text(txt)
#if tspan.node().getComputedTextLength() > width)
#alert(tspan.node().getComputedTextLength())
#fill = d3.scale.category10()
#nodes = d3.range(100).map(Object)
#
#vis = d3.select("#chart").append("svg:svg")
#.attr("width", w)
#.attr("height", h)
#
#force = d3.layout.force()
#.nodes(nodes)
#.links([])
#.size([w, h])
#.start()
#
#node = vis.selectAll("circle.node")
#.data(nodes)
#.enter().append("svg:circle")
#.attr("class", "node")
#.attr("cx", (d) -> d.x)
#.attr("cy", (d) -> return d.y)
#.attr("r", 8)
#.style("fill", (d, i) -> fill(i & 3) )
#.style("stroke", (d, i) -> d3.rgb(fill(i & 3)).darker(2) )
#.style("stroke-width", 1.5)
#.call(force.drag)
#
#vis.style("opacity", 1e-6)
#.transition()
#.duration(1000)
#.style("opacity", 1)
#
#force.on "tick", (e) ->
## Push different nodes in different directions for clustering.
#  k = 6 * e.alpha
#  nodes.forEach (o, i) ->
#    o.x += if i & 2 then k else -k
#    o.y += if i & 1 then k else -k
#  node.attr("cx", (d) -> d.x )
#  .attr("cy", (d) -> d.y )
#
#
#d3.select("body").on "click", () ->
#  nodes.forEach (o, i) ->
#    o.x += (Math.random() - 0.5) * 40
#    o.y += (Math.random() - 0.5) * 40
#  force.resume()