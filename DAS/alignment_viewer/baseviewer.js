// Generated by CoffeeScript 1.12.2
(function() {
  var contig_height, fill, h, primer_fill, query_height, spacer, svg, tooltip, w, xpadding, y;

  w = 960;

  h = 5000;

  svg = d3.select('body').append('svg').attr('width', w).attr('height', h);

  xpadding = 50.0;

  y = 50;

  spacer = 12;

  query_height = 10;

  contig_height = 7;

  tooltip = d3.select("body").append("div").attr("class", "tooltip");

  d3.json('data.json', function(data) {
    var query_len, xScale, yScale, yend;
    query_len = data.contigs[0].query._Region__length;
    xScale = d3.scaleLinear().domain([0, query_len]).range([0 + xpadding, w - xpadding]);
    yScale = d3.scaleLinear().domain([0, data.contigs.length]).range([y, data.contigs.length * spacer]);
    yend = yScale(data.contigs.length + 2);
    svg.append('text').attr('x', xScale(query_len / 2.0)).attr('y', yScale(-2)).attr("text-anchor", "middle").text(function(d) {
      return data.meta.query + " (" + data.meta.query_length + " bp)";
    });
    svg.append('rect').attr('x', xScale(0)).attr('width', xScale(query_len)).attr('y', yScale(0)).attr('height', query_height).attr('opacity', 0.5).attr('fill', 'blue');
    svg.append('line').attr('x1', xScale(query_len / 2.0)).attr('y1', yScale(0)).attr('x2', xScale(query_len / 2.0)).attr('y2', yend).style("stroke-dasharray", "3, 3").style("stroke", 'black').style("stroke-width", 1.5);
    svg.selectAll('rects').data(data.contigs).enter().append('rect').attr('x', function(d) {
      return xScale(d.query._Region__start);
    }).attr('y', function(d, i) {
      return yScale(i + 2);
    }).attr('width', function(d) {
      return xScale(d.query._Region__end) - xScale(d.query._Region__start);
    }).attr('height', contig_height).attr('fill', fill).attr('opacity', 0.9).on("mouseover", function(d) {
      var coordinates;
      coordinates = d3.mouse(this);
      d3.select(this).style("fill", 'red');
      return tooltip.html('<b>Subject: </b>' + d.subject.name + '<br>' + '<b>QRange:\t</b> ' + d.query._Region__start + '-' + d.query._Region__end + '<br>' + '<b>Type:\t</b>' + d.contig_type + '<br>' + '<b>Id:\t</b>' + d.contig_id + '<br>' + '<b>SRange:\t</b>' + d.subject._Region__start + '-' + d.subject._Region__end + '<br>' + '<b>Ends:\t</b>' + d.start_label + ' ' + d.end_label).style("visibility", "visible");
    }).on("mousemove", function() {
      return tooltip.style("top", (d3.event.pageY - 10) + "px").style("left", (d3.event.pageX + 10) + "px");
    }).on("mouseout", function(d) {
      d3.select(this).style('fill', fill);
      return tooltip.style("visibility", "hidden");
    });
    svg.selectAll('forward_primers').data(data.contigs).enter().append('rect').attr('x', function(d) {
      return xScale(d.query._Region__start);
    }).attr('y', function(d, i) {
      return yScale(i + 2);
    }).attr('width', 2).attr('height', contig_height).attr('fill', function(d) {
      var c;
      c = 'black';
      if (d.start_label === 'new_primer') {
        c = 'red';
      }
      return c;
    });
    svg.selectAll('forward_primers').data(data.contigs).enter().append('rect').attr('x', function(d) {
      return xScale(d.query._Region__end - 2);
    }).attr('y', function(d, i) {
      return yScale(i + 2);
    }).attr('width', 2).attr('height', contig_height).attr('fill', function(d) {
      var c;
      c = 'black';
      if (d.end_label === 'new_primer') {
        c = 'red';
      }
      return c;
    });
    d3.json('primer_data.json', function(primerdata) {
      return svg.selectAll('primer_rects').data(primerdata.contigs).enter().append('rect').attr('x', function(d) {
        return xScale(d.query._Region__start);
      }).attr('y', function(d, i) {
        return yScale(0);
      }).attr('width', function(d) {
        return xScale(d.query._Region__end) - xScale(d.query._Region__start);
      }).attr('height', query_height * 0.75).attr('fill', "red").attr('opacity', 1.35);
    });
    return d3.json('primer_data.json', function(primerdata) {
      return svg.selectAll('primer_text').data(primerdata.contigs).enter().append('text').append('text').attr('x', function(d) {
        return xScale(d.query._Region__start);
      }).attr('y', function(d, i) {
        return yScale(-2);
      }).text('OK');
    });
  });

  primer_fill = function(d) {
    var c;
    c = 'black';
    if (d.start_label === 'new_primer') {
      c = 'red';
    }
    return c;
  };

  fill = function(d) {
    if (d.contig_type === 'blast') {
      return 'purple';
    }
    return 'orange';
  };

}).call(this);

//# sourceMappingURL=baseviewer.js.map
