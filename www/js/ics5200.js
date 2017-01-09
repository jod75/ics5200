$(document).ready(function () {
	var table = $('#similarityTable').DataTable();
 
    $.ajax({
        url: "http://hadoop1:5432/getTestMolRegNos",
        type: "get",
        datatype: "json",
        success: function (response) {
            mols = JSON.parse(response)
            // fill drop down
            var dropdown = $("#testMolecules");
            dropdown.empty(); 
            for (var i = 0; i < mols.length; i++) {
                dropdown.append("<option>" + mols[i] + "</option>");
            }          
        }
    });
 
	$('#similarityTable tbody').on( 'click', 'tr', function () {
        if (table.row(this).length > 0) {
            console.log( 'Row index: '+table.row( this ).index() );
            $.ajax({
                url: "http://hadoop1:5432/smilesToSVG/" + table.row(this).data()[2],
                type: "get",
                datatype: "json",
                success: function(response) {                    
                    document.getElementById('molecule').innerHTML = response                    
                }
            });
        }
	} );
		
	$("#userQueryButton").click(function () {		
        $.ajax({
            url: "http://hadoop1:5432/testSmilesToSVG/" + $("#testMolecules").val(),
            type: "get",
            datatype: "json",
            success: function(response) {                    
                document.getElementById('queryMolecule').innerHTML = response       
                document.getElementById('molecule').innerHTML = ""              
            }
        });
		$.ajax({
			url: "http://hadoop1:5432/1/sample/" + $("#testMolecules").val(),
			type: "get",
			datatype: "json",			
			success: function (response) {				
				res = JSON.parse(response);
                table.rows().clear();
				table.rows.add(res).draw();
			}
		});
	});
});

