
var ridge_plot = {
    viewBox:{ left:0, top: 0, width:500, height:400 },
    margin:{top: 10, right: 10, bottom: 30, left: 80}
}

function setup_ridge_plot(data, container) {

    // Setup width and height
    ridge_plot.width = ridge_plot.viewBox.width -
        ridge_plot.margin.right - ridge_plot.margin.left;
    ridge_plot.height = ridge_plot.viewBox.height -
        ridge_plot.margin.top - ridge_plot.margin.bottom;

    ridge_plot.svg = container.attr('viewBox', Object.values(ridge_plot.viewBox))
                              .append('g')
                              .attr("transform", `translate(${ridge_plot.margin.left} ${ridge_plot.margin.top})`)


    // Setup scales for plot
    ridge_plot.scale = {}
    ridge_plot.scale.band = d3.scaleBand()
                              .range([0, ridge_plot.height])
                              .domain(data.map(d => d.name))
                              .padding(0.4);
    ridge_plot.scale.x = d3.scaleLinear()
                           .domain([0, data[0].good.length-2])
                           .range([0, ridge_plot.width])
                           .nice()
    ridge_plot.scale.y = d3.scaleLinear()
                           .domain([0, d3.max(data.map((d)=>d3.mean(d.good)*2))])
                           .range([ridge_plot.scale.band.bandwidth(),0])


    // add x axis
    ridge_plot.svg.append("g")
              .attr("transform", `translate(0,${ridge_plot.height})`)
              .call(d3.axisBottom(ridge_plot.scale.x)
                      .tickFormat( (d, i) => d3.format(".1s")(d)))

    // setup draw functions for output and area
    ridge_plot.draw = { outline:d3.line()
                                  .curve(d3.curveBasis)
                                  .x((d,i) => ridge_plot.scale.x(i))
                                  .y(d => ridge_plot.scale.y(d)),
                        area:d3.area()
                               .curve(d3.curveBasis)
                                  .y1(d => ridge_plot.scale.y(d))
                                  .y0(ridge_plot.scale.y(0))
                                  .x((d,i) => ridge_plot.scale.x(i))
    };


    // add disto plot groups
    ridge_plot.distros = ridge_plot.svg
                                   .selectAll(".distro")
                                   .data(data)
                                   .enter()
                                   .append("g")
                                   .classed("distro", true)
                                   .attr("transform",
                                         d=>`translate(0, ${ridge_plot.scale.band(d.name)})`)

    // add distro outlines
    ridge_plot.distros.append("path")
              .classed("line", true)
              .attr("d", (d)=>ridge_plot.draw.outline(d.good.slice(1,-1)))
    // add distro area
    ridge_plot.distros.append("path")
              .classed("area", true)
              .attr("d", (d)=>ridge_plot.draw.area(d.good.slice(1,-1)))
    // add distro text
    ridge_plot.distros.append("text")
           .text(d => d.name)
           .attr('y', ridge_plot.scale.band.bandwidth)
           .attr('x', -ridge_plot.margin.left/2)


    // setup update function
    ridge_plot.update_counts = function (){
        const selection = d3.event.selection;

      ridge_plot.extent = selection
            .map(ridge_plot.scale.x.invert)
            .map(Math.round)



        // Get total number of fragments under selection
        sizes = data.map(function(d) {
            return { name: d.name ,
                     total:d3.sum(d.good.slice(ridge_plot.extent[0], ridge_plot.extent[1] + 1)) }
        })

        // Sort plot bands by selection total
        sizes.sort((a,b) => b.total - a.total)
        ridge_plot.scale.band.domain(sizes.map(d=>d.name))


        // update labels with new totals
        ridge_plot.svg.selectAll(".size_label")
                  .data(sizes)
                  .join("text")
                  .classed("size_label", true)
                  .text(d => d.total )
                  .attr("y", d => ridge_plot.scale.band(d.name) )
                  .attr("x", d => d3.mean(selection))

        // update distro plots to new sort
        ridge_plot.distros
                  .transition()
                  .attr("transform", d=>`translate(0, ${ridge_plot.scale.band(d.name)})`)

      /* update table if the code has loaded */
      if(table_update) table_update();

        if (!d3.event.sourceEvent || !selection) return;
      d3.select(this).transition().call(ridge_plot.brush.func.move, ridge_plot.extent.map(ridge_plot.scale.x))

    }

    // Update labels
    ridge_plot.update_labels = function(){
        const selection = d3.event.selection;

      ridge_plot.extent = selection
            .map(ridge_plot.scale.x.invert)
            .map(Math.round)

        ridge_plot.brush.labels.left
                  .attr('x', ridge_plot.scale.x(ridge_plot.extent[0]))
                  .text(ridge_plot.extent[0])
        ridge_plot.brush.labels.right
                  .attr('x', ridge_plot.scale.x(ridge_plot.extent[1]))
                  .text(ridge_plot.extent[1])

        if (!d3.event.sourceEvent || !selection) return;
    }


    ridge_plot.brush = {func:d3.brushX()
                         .extent([[0, 0], [ridge_plot.width, ridge_plot.height]])
                         .on("end", ridge_plot.update_counts)
                         .on("brush", ridge_plot.update_labels)}


    ridge_plot.brush.element = ridge_plot.svg.append("g")
                                  .call(ridge_plot.brush.func)
                                  .call(g => g.select(".overlay"))


    ridge_plot.brush.labels = {
     left: ridge_plot.svg.append('text')
              .attr('id', 'labelleft')
              .attr('x', 0)
              .attr('y', 0),
     right: ridge_plot.svg.append('text')
              .attr('id', 'labelright')
              .attr('x', 0)
              .attr('y', 0)
 }


    ridge_plot.brush.element.call(ridge_plot.brush.func.move,
                                  [150, 700].map(ridge_plot.scale.x))

}
