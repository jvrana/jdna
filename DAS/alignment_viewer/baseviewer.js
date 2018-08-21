// Generated by CoffeeScript 1.12.2
(function() {
  var fill, h, svg, tooltip, w, xpadding;

  w = 960;

  h = 5000;

  svg = d3.select('body').append('svg').attr('width', w).attr('height', h);

  xpadding = 50.0;

  tooltip = d3.select("body").append("div").attr("class", "tooltip");

  d3.json('data.json', function(data) {
    var query_len, xScale, yScale;
    query_len = data.contigs[0].query_length;
    xScale = d3.scaleLinear().domain([0, query_len]).range([0 + xpadding, w - xpadding]);
    yScale = d3.scaleLinear().domain([0, data.contigs.length]).range([120, 500]);
    svg.append('rect').attr('x', xScale(0)).attr('width', xScale(query_len)).attr('y', 100).attr('height', 10).attr('opacity', 0.5).attr('fill', 'blue');
    svg.selectAll('rects').data(data.contigs).enter().append('rect').attr('x', function(d) {
      return xScale(d.q_start);
    }).attr('y', function(d, i) {
      return 120 + i * 10;
    }).attr('width', function(d) {
      return xScale(d.q_end) - xScale(d.q_start);
    }).attr('height', 7).attr('fill', fill).attr('opacity', 0.9).on("mouseover", function(d) {
      var coordinates;
      coordinates = d3.mouse(this);
      d3.select(this).style("fill", 'red');
      return tooltip.text(d.subject_acc + '<br>' + d.q_start + ', ' + d.q_end + '<br>' + d.contig_type + ' ' + d.contig_id).style("visibility", "visible");
    }).on("mousemove", function() {
      return tooltip.style("top", (d3.event.pageY - 10) + "px").style("left", (d3.event.pageX + 10) + "px");
    }).on("mouseout", function(d) {
      d3.select(this).style('fill', fill);
      return tooltip.style("visibility", "hidden");
    });
    svg.selectAll('rects').data(data.contigs).enter().append('rect').attr('x', function(d) {
      return xScale(d.q_start);
    }).attr('y', function(d, i) {
      return 100;
    }).attr('width', 1).attr('height', h).attr('fill', 'black').attr('opacity', 0.1);
    svg.selectAll('rects').data(data.contigs).enter().append('rect').attr('x', function(d) {
      return xScale(d.q_end);
    }).attr('y', function(d, i) {
      return 100;
    }).attr('width', 1).attr('height', h).attr('fill', 'blue').attr('opacity', 0.1);
    return d3.json('primer_data.json', function(primerdata) {
      return svg.selectAll('primer_rects').data(primerdata.contigs).enter().append('rect').attr('x', function(d) {
        return xScale(d.q_start);
      }).attr('y', function(d, i) {
        return 100;
      }).attr('width', 1).attr('height', h).attr('fill', "red").attr('opacity', 0.35);
    });
  });

  fill = function(d) {
    if (d.contig_type === 'contig') {
      return 'black';
    }
    if (d.contig_type === 'product') {
      return 'purple';
    }
    if (d.contig_type === 'gap') {
      return 'blue';
    }
    return 'orange';
  };

}).call(this);

//# sourceMappingURL=baseviewer.js.map