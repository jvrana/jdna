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
y = 50
spacer = 12
query_height = 10
contig_height = 7

tooltip = d3.select("body").append("div")
  .attr("class", "tooltip")


d3.json('data.json', (data) ->



  query_len = data.contigs[0].query_length

  xScale = d3.scaleLinear()
    .domain([0, query_len])
    .range([0+xpadding, w-xpadding])

  yScale = d3.scaleLinear()
    .domain([0, data.contigs.length])
    .range([y, data.contigs.length * spacer])

  yend = yScale(data.contigs.length + 2)

  svg.append('text')
    .attr('x', xScale(query_len/2.0))
    .attr('y', yScale(-2))
    .attr("text-anchor", "middle")
    .text( (d) -> data.meta.query + " (" + data.meta.query_length + " bp)")

  svg.append('rect')
    .attr('x', xScale(0))
    .attr('width', xScale(query_len))
    .attr('y', yScale(0))
    .attr('height', query_height)
    .attr('opacity', 0.5)
    .attr('fill', 'blue')

  svg.append('line')
    .attr('x1', xScale(query_len/2.0))
    .attr('y1', yScale(0))
    .attr('x2', xScale(query_len/2.0))
    .attr('y2', yend)
    .style("stroke-dasharray", ("3, 3"))
    .style("stroke", 'black' )
    .style("stroke-width", 1.5)



  svg.selectAll('rects')
    .data(data.contigs)
    .enter()
    .append('rect')
    .attr('x', (d) -> xScale(d.q_start))
    .attr('y', (d,i) -> yScale(i+2))
    .attr('width', (d) -> xScale(d.q_end) - xScale(d.q_start))
    .attr('height', contig_height)
    .attr('fill', fill)
    .attr('opacity', 0.9)
    .on("mouseover", (d) ->
      coordinates = d3.mouse(this)
      d3.select(this).style("fill", 'red')
      tooltip
        .html('<b>Subject: </b>' + d.subject_acc + '<br>' +
          '<b>Range:\t</b> ' + d.q_start + '-' + d.q_end + '<br>' +
          '<b>Type:\t</b>' + d.contig_type + '<br>' +
          '<b>Id:\t</b>' + d.contig_id)
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
    .attr('y', (d,i) -> yScale(1))
    .attr('width', 1)
    .attr('height', yend - yScale(1))
    .attr('fill', 'black')
    .attr('opacity', 0.1)

  svg.selectAll('rects')
    .data(data.contigs)
    .enter()
    .append('rect')
    .attr('x', (d) -> xScale(d.q_end))
    .attr('y', (d,i) -> yScale(1))
    .attr('width', 1)
    .attr('height', yend - yScale(1))
    .attr('fill', 'blue')
    .attr('opacity', 0.1)

  d3.json('primer_data.json', (primerdata) ->
      svg.selectAll('primer_rects')
        .data(primerdata.contigs)
        .enter()
        .append('rect')
        .attr('x', (d) -> xScale(d.q_start))
        .attr('y', (d,i) -> yScale(0))
        .attr('width', (d) -> xScale(d.q_end) - xScale(d.q_start))
        .attr('height', query_height * 0.75)
        .attr('fill', "red")
        .attr('opacity', 1.35)
  )

  d3.json('primer_data.json', (primerdata) ->
      svg.selectAll('primer_text')
        .data(primerdata.contigs)
        .enter()
        .append('text')
        .append('text')
        .attr('x', (d) -> xScale(d.q_start))
        .attr('y', (d,i) -> yScale(-2))
        .text('OK')

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