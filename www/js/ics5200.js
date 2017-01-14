function sortNumber(a,b) {
    return a - b;
}

$(document).ready(function () {
	var molSimTable = $('#ligandSimilarityTable').DataTable();
    var protSimTable = $('#proteinSimilarityTable').DataTable();
    var testLigandDataTable = $('#testLigandDataTable').DataTable();
    var testProteinDataTable = $('#testProteinDataTable').DataTable();
    
    $("#tabs").tabs();    
    
    // populate ligand test data table
    $.ajax({
        url: "http://hadoop1:5432/getTestLigandsDS",
        type: "get",
        datatype: "json",
        success: function (response) {
            var res = JSON.parse(response);
            testLigandDataTable.rows().clear();   
            testLigandDataTable.rows.add(res).draw();  
        }
    });
    
    // populate protein test data table
    $.ajax({
        url: "http://hadoop1:5432/getTestProteinsDS",
        type: "get",
        datatype: "json",
        success: function (response) {
            var res = JSON.parse(response);
            testProteinDataTable.rows().clear();                           
            testProteinDataTable.rows.add(res).draw();  
        }
    });
     
	$('#ligandSimilarityTable tbody').on('click', 'tr', function () {
        if (molSimTable.row(this).length > 0) {            
            document.getElementById('knownLigandMolRegNo').innerText = molSimTable.row(this).data()[0]; 
            $.ajax({
                url: "http://hadoop1:5432/getSmilesSVG/" + molSimTable.row(this).data()[0],
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
            var molRegNo = testLigandDataTable.row(this).data()[6];
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
                var targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    testLigandBindings.innerHTML += ("<li><samp>" + targets[i] + "</samp></li>");
                }   
                testLigandBindings.innerHTML += ("</ul>");
            }
            });
        }
	} );

    $('#testProteinDataTable tbody').on('click', 'tr', function () {
        if (testProteinDataTable.row(this).length > 0) { 
            var compId = testProteinDataTable.row(this).data()[0];
            document.getElementById("selectedTestProtein").innerText = compId;

            /*$.ajax({
            url: "http://hadoop1:5432/getProteinTestBindingLigands/" + compId,
            type: "get",
            datatype: "json",
            success: function(response) {  
                var testLigandBindings = document.getElementById("testProteinBindings")              
                testLigandBindings.innerHTML = "<ul>"
                targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    testLigandBindings.innerHTML += ("<li><samp>" + targets[i] + "</samp></li>");
                }   
                testLigandBindings.innerHTML += ("</ul>");
            }
            });*/
        }
	} );
    
    $('#testLigandDataTable tbody').on('dblclick', 'tr', function () {
        if (testLigandDataTable.row(this).length > 0) {   
            var molRegNo = testLigandDataTable.row(this).data()[6]
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
                queryLigandMolRegNo.innerHTML = molRegNo + " ( "
                var targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    queryLigandMolRegNo.innerHTML += (targets[i][0] + " ");
                }   
                queryLigandMolRegNo.innerHTML += ")";
            }
            });

            $.ajax({
			url: "http://hadoop1:5432/doLigandExperiment/" + molRegNo,
			type: "get",
			datatype: "json",			
			success: function (response) {				
				var res = JSON.parse(response);
                molSimTable.rows().clear();
                molSimTable.column(2).visible(false);
                molSimTable.column(3).visible(false);
                molSimTable.column(4).visible(false);
				molSimTable.rows.add(res).draw();

                var knownLigandsUniqueCompId = document.getElementById("knownLigandsUniqueCompId");
                knownLigandsUniqueCompId.innerHTML = "( ";
                var uniqueCompIds = molSimTable.column(10).data().unique().toArray().sort(sortNumber);
                for (var i = 0; i < uniqueCompIds.length; i++) {
                    knownLigandsUniqueCompId.innerHTML += (uniqueCompIds[i] + " ");
                }
                knownLigandsUniqueCompId.innerHTML += ")";
                $('#tabs').tabs({ active: 2 });
			}
		    });
        }
	} );
    
    /*$(document).ajaxSend(function(event, jqxhr, settings){
        if (settings.url.includes("/doLigandExperiment/") || settings.url.includes("/getLigandTestTargets/")) {
            $.LoadingOverlay("show", {image: "./images/molecules-bonding-animation-5.gif", size: "100%", maxSize: 240, resizeInterval: 0});
        }        
    });

    $(document).ajaxComplete(function(event, jqxhr, settings){ 
        if (settings.url.includes("/doLigandExperiment/" ) || settings.url.includes("/getLigandTestTargets/") ) {
            $.LoadingOverlay("hide");
        }        
    });*/
});

