$(document).ready(function () {
	var molSimTable = $('#moleculeSimilarityTable').DataTable();
    var protSimTable = $('#proteinSimilarityTable').DataTable();
    var testDataTable = $('#testdataTable').DataTable();
    
    $("#tabs").tabs();
 
    // fill dropdown list
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
    
    // populate test data table
    $.ajax({
        url: "http://hadoop1:5432/getTestDataset",
        type: "get",
        datatype: "json",
        success: function (response) {
            res = JSON.parse(response);
            testDataTable.rows().clear();   
            testDataTable.column(0).visible(false);            
            testDataTable.rows.add(res).draw();  
        }
    });
    
    // binds to tab events
    $('#tabs').bind('tabsactivate', function(event, ui) { 
      switch (ui.newPanel[0].id){
        case "tabs-1": 
            console.log("0")
            break;
        case "tabs-2": 
            console.log("1")
            break;
        case "tabs-3": 
            console.log("2")
            break;
      }
    });
 
	$('#moleculeSimilarityTable tbody').on('click', 'tr', function () {
        if (molSimTable.row(this).length > 0) {
            console.log( 'Row index: ' + molSimTable.row(this).index() );
            document.getElementById('selectedExperimentMolecule').innerText = molSimTable.row(this).data()[2]; 
            $.ajax({
                url: "http://hadoop1:5432/smilesToSVG/" + molSimTable.row(this).data()[2],
                type: "get",
                datatype: "json",
                success: function(response) {                    
                    document.getElementById('molecule').innerHTML = response                    
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
    
    $('#testdataTable tbody').on('click', 'tr', function () {
        if (testDataTable.row(this).length > 0) { 
            document.getElementById("selectedTestMolecule").innerText = testDataTable.row(this).data()[2];
            $.ajax({
                url: "http://hadoop1:5432/testSmilesToSVG/" + testDataTable.row(this).data()[2],
                type: "get",
                datatype: "json",
                success: function(response) {                    
                    document.getElementById('testmolecule').innerHTML = response                    
                }
            });
        }
	} );
    
    $('#testdataTable tbody').on('dblclick', 'tr', function () {
        if (testDataTable.row(this).length > 0) {   
            var molreg = testDataTable.row(this).data()[2]
            document.getElementById('queryProtein').innerText = testDataTable.row(this).data()[8]
            var dd = document.getElementById('testMolecules');
            for (var i = 0; i < dd.options.length; i++) {
                if (dd.options[i].text === molreg) {
                    dd.selectedIndex = i;
                    break;
                }
            }            
            $("#userQueryButton").click();
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
            url: "http://hadoop1:5432/getTestTargets/" + $("#testMolecules").val(),
            type: "get",
            datatype: "json",
            success: function(response) {                    
                document.getElementById('targets').innerHTML = "<ul>"
                targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    document.getElementById('targets').innerHTML += ("<li>" + targets[i] + "</li>");
                }   
                document.getElementById('targets').innerHTML += ("</ul>");
            }
        });
		$.ajax({
			url: "http://hadoop1:5432/doExperiment/" + $("#testMolecules").val(),
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

