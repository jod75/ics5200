$(document).ready(function () {
	var molSimTable = $('#ligandSimilarityTable').DataTable();
    var protSimTable = $('#proteinSimilarityTable').DataTable();
    var testLigandDataTable = $('#testLigandDataTable').DataTable();
    var testProteinDataTable = $('#testProteinDataTable').DataTable();
    
    $("#tabs").tabs();    
    
    // populate test data table
    $.ajax({
        url: "http://hadoop1:5432/getTestLigandsDS",
        type: "get",
        datatype: "json",
        success: function (response) {
            res = JSON.parse(response);
            testLigandDataTable.rows().clear();   
            testLigandDataTable.column(0).visible(false);            
            testLigandDataTable.rows.add(res).draw();  
        }
    });
    
    $.ajax({
        url: "http://hadoop1:5432/getTestProteinsDS",
        type: "get",
        datatype: "json",
        success: function (response) {
            res = JSON.parse(response);
            testProteinDataTable.rows().clear();   
            testProteinDataTable.column(0).visible(false);            
            testProteinDataTable.rows.add(res).draw();  
        }
    });
     
	$('#ligandSimilarityTable tbody').on('click', 'tr', function () {
        if (molSimTable.row(this).length > 0) {
            console.log( 'Row index: ' + molSimTable.row(this).index() );
            document.getElementById('knownLigandMolRegNo').innerText = molSimTable.row(this).data()[2]; 
            $.ajax({
                url: "http://hadoop1:5432/getSmilesSVG/" + molSimTable.row(this).data()[2],
                type: "get",
                datatype: "json",
                success: function(response) {                    
                    document.getElementById('knownLigandSVG').innerHTML = response                    
                }
            });
        }
	} );
    
    $('#proteinSimilarityTable tbody').on('click', 'tr', function () {
        if (protSimTable.row(this).length > 0) {
            console.log( 'Row index: ' + protSimTable.row(this).index() );
            document.getElementById('proteinAlignment').innerText =  protSimTable.row(this).data()[13].replace(/ /g, '\u00a0') 
        }
	} );
    
    $('#testLigandDataTable tbody').on('click', 'tr', function () {
        if (testLigandDataTable.row(this).length > 0) { 
            var molRegNo = testLigandDataTable.row(this).data()[2];
            document.getElementById("testLigandMolRegNo").innerText = molRegNo;
            $.ajax({
                url: "http://hadoop1:5432/getSmilesSVG/" + molRegNo,
                type: "get",
                datatype: "json",
                success: function(response) {                    
                    document.getElementById('testLigandSVG').innerHTML = response                    
                }
            });

            $.ajax({
            url: "http://hadoop1:5432/getLigandTestTargets/" + molRegNo,
            type: "get",
            datatype: "json",
            success: function(response) {  
                var testLigandBindings = document.getElementById("testLigandBindings")              
                testLigandBindings.innerHTML = "<ul>"
                targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    testLigandBindings.innerHTML += ("<li><samp>" + targets[i] + "</samp></li>");
                }   
                testLigandBindings.innerHTML += ("</ul>");
            }
            });
        }
	} );
    
    $('#testLigandDataTable tbody').on('dblclick', 'tr', function () {
        if (testLigandDataTable.row(this).length > 0) {   
            var molRegNo = testLigandDataTable.row(this).data()[2]
            document.getElementById("queryLigandMolRegNo").innerText = molRegNo;
            $.ajax({
                url: "http://hadoop1:5432/getSmilesSVG/" + molRegNo,
                type: "get",
                datatype: "json",
                success: function(response) {                    
                    document.getElementById('queryLigandSVG').innerHTML = response                    
                }
            });

            $.ajax({
            url: "http://hadoop1:5432/getLigandTestTargets/" + molRegNo,
            type: "get",
            datatype: "json",
            success: function(response) {  
                var queryLigandMolRegNo = document.getElementById("queryLigandMolRegNo")              
                queryLigandMolRegNo.innerText = molRegNo + " ( "
                targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    queryLigandMolRegNo.innerText += (targets[i] + " ");
                }   
                queryLigandMolRegNo.innerText += ")";
            }
            });

            $.ajax({
			url: "http://hadoop1:5432/doLigandExperiment/" + molRegNo,
			type: "get",
			datatype: "json",			
			success: function (response) {				
				res = JSON.parse(response);
                molSimTable.rows().clear();
                molSimTable.column(0).visible(false);
                molSimTable.column(9).visible(false);
				molSimTable.rows.add(res).draw();
                $('#tabs').tabs({ active: 2 });
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
			url: "http://hadoop1:5432/doProteinExperiment/" + document.getElementById('queryProtein').innerText,
			type: "get",
			datatype: "json",			
			success: function (response) {				
				res = JSON.parse(response);
                protSimTable.rows().clear();
                protSimTable.column(0).visible(false);
                protSimTable.column(12).visible(false);
                protSimTable.column(13).visible(false);
				protSimTable.rows.add(res).draw();
                $('#tabs').tabs({ active: 3 });
			}
		});
	});
    
    $(document).ajaxSend(function(event, jqxhr, settings){
        if (settings.url.includes("/doExperiment/" )) {
            $.LoadingOverlay("show", {image: "./images/molecules-bonding-animation-5.gif", size: "100%", maxSize: 240, resizeInterval: 0});
        }        
    });

    $(document).ajaxComplete(function(event, jqxhr, settings){
        if (settings.url.includes("/doExperiment/" )) {
            $.LoadingOverlay("hide");
        }        
    });
});

