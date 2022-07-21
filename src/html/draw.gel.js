
var gel_plot = {
    viewBox:{ left:0, top: 0, width:460, height:400 },
    margin:{top: 100, right: 10, bottom: 10, left: 30}
}

function correct_enzyme_name (d) {
  var el = d3.select(this);
  var words=d3.select(this).text().split(/(\b[A-Z][a-z]*)/);

  el.text('');

  for (var i = 1; i < words.length; i++) {
    var tspan = el.append('tspan').text(words[i]);
    if((i%2) == 1)
      tspan.attr("style", "font-style:italic");
  }
 }

function setup_gel_plot(data, container) {

    // Setup width and height
    gel_plot.width = gel_plot.viewBox.width -
        gel_plot.margin.right - gel_plot.margin.left;
    gel_plot.height = gel_plot.viewBox.height -
        gel_plot.margin.top - gel_plot.margin.bottom;

    gel_plot.svg = container.attr('viewBox', Object.values(gel_plot.viewBox))
                            .append('g')
                            .attr("transform", `translate(${gel_plot.margin.left} ${gel_plot.margin.top})`)

    gel_plot.svg
            .append('rect')
            .attr('width', '100%')
            .attr('height', '100%')

    // Setup scales for plot
    gel_plot.scale = {}
    gel_plot.scale.x = d3.scaleBand()
                         .range([0, gel_plot.width])
                         .domain(data.map(d => d.name))
                         .padding(0.1);
    gel_plot.scale.y = d3.scaleLog()
                           .base(Math.E)
                           .domain([100, data[0].all.length])
                           .range([gel_plot.height, 0])
                           .clamp(true);
    gel_plot.scale.alpha = d3.scaleSqrt([0,data[0].genome_size], [0,0.70])

    // add y axis
    gel_plot.svg
            .append("g")
            .call(d3.axisLeft(gel_plot.scale.y)
                    .tickValues([100, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000,10000])
                    .tickFormat( (d, i) => d3.format(".2s")(d)))

   // add x axis
    gel_plot.svg.append("g")
              .call(d3.axisTop(gel_plot.scale.x))
      .selectAll("text")	
        .style("text-anchor", "start")
        .attr("dx", ".8em")
        .attr("dy", ".15em")
        .attr("transform", "rotate(-65)")
        .each(correct_enzyme_name)
    // add lane plot groups
    gel_plot.lanes = gel_plot.svg
                             .selectAll(".lanes")
                             .data(data)
                             .enter()
                             .append("g")
                             .classed("distro", true)
                             .attr("transform",
                                   d=>`translate(${gel_plot.scale.x(d.name)},0)`)
    // add data for gel lanes
    gel_plot.lanes.selectAll("rect")
            .data((d) => d.all)
            .enter()
            .append('rect')
            .classed("entry", true)
            .attr('y', (d,i) => gel_plot.scale.y(i)-2.5)
            .attr('width', gel_plot.scale.x.bandwidth)
            .attr('opacity', (d,i) => gel_plot.scale.alpha(d*i))
            .attr('height', 5)
            .attr('fill', 'white');
    

}
