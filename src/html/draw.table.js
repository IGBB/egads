var tbody;


function setup_table(data, table) {

    tbody = table.append('tbody')

    table_update();
}

var table_update = function () {

    if(!tbody) return;


    var row_data = data.map( function (row) {
        total_selected = d3.sum(row.good.slice(ridge_plot.extent[0],
                                               ridge_plot.extent[1] + 1));

        methyl_frag = Math.round(document.getElementById('methyl').value/100 *
                                 total_selected);

        mutat_frag = row.mutation;

        miss_frag = Math.round(document.getElementById('missing').value/100 *
                               total_selected);

        read_depth = document.getElementById('depth').value *
            ( total_selected -
              methyl_frag -
              mutat_frag -
              miss_frag );

        samp_per_lane = Math.round(document.getElementById('lane').value *
                                   1000000 /read_depth);

        reduction =  d3.sum(row.good.map((d,i) => d*i)
                              .slice(ridge_plot.extent[0],
                                     ridge_plot.extent[1] + 1))/ row.genome_size

      if(read_depth < 0){
        read_depth = "-";
        samp_per_lane = "-";
        reduction = "-"
      }




        var compat = ""
        if(row.compat == 0) compat = " &#9940;"
        if(row.compat == 1) compat = " &#9989;"
        if(row.compat == 2) compat = " &#10067;"

        return [compat, // enzyme compatability
                row.name.replace(/\b[A-Z][a-z]*/g, "<i>$&</i>"), // enzyme pair
                d3.sum(row.good), // total number of good fragments
                total_selected, // total number of good fragments inside selected size
                methyl_frag, // methylation reduction
                mutat_frag, // mutation reduction
                miss_frag, // allowed missing
                read_depth,  // total read needed
                samp_per_lane, // Number of samples per lane
                reduction.toLocaleString("US", {style: "percent"}) // reduction rate
               ];
    })


    row_data.sort((a,b) => b[3] - a[3]);


    // create a row for each object in the data
    var rows = tbody.selectAll("tr")
                    .data(row_data)
                    .join("tr");

    // create a cell in each row for each column
    var cells = rows.selectAll("td")
                    .data((d) => d)
                    .join("td")
                    .html((d,i) => (i > 0 && i < 8)? d.toLocaleString() : d);

}
