function sortNumber(a, b) {
    return a - b;
}

function fillDropDown(dropDownElement, values, defaultSelectionIndex) {
    var sel = document.getElementById(dropDownElement);
    var fragment = document.createDocumentFragment();

    values.forEach(function(valueItem, index) {
        var opt = document.createElement("option");
        opt.innerHTML = valueItem;
        opt.value = valueItem;
        fragment.appendChild(opt);
    });

    sel.appendChild(fragment);
    sel.selectedIndex = defaultSelectionIndex;
}

function setTextBoxText(textBoxElement, value){
    var sel = document.getElementById(textBoxElement);
    sel.value = value;
}

$(document).ready(function () {
	var molSimTable = $('#ligandSimilarityTable').DataTable();
    var protSimTable = $('#proteinSimilarityTable').DataTable();
    var testLigandDataTable = $('#testLigandDataTable').DataTable();
    var testProteinDataTable = $('#testProteinDataTable').DataTable();
    var ligSimilarityAlgorithms = ["None", "Tanimoto", "Dice"];
    var ligSimilarityFingerprints = ["Morgan", "MACCS"];
    
    $("#tabs").tabs();    

    fillDropDown("ligSim1", ligSimilarityAlgorithms, 1);
    fillDropDown("ligSim2", ligSimilarityAlgorithms, 0);
    fillDropDown("ligFP1", ligSimilarityFingerprints, 0);
    fillDropDown("ligFP2", ligSimilarityFingerprints, 0);
    setTextBoxText("ligTH1", "0.8");
    setTextBoxText("ligTH2", "0.5");
    
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
                url: "http://hadoop1:5432/getSmilesSVG/" + molSimTable.row(this).data()[0] + "/mol",
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
            document.getElementById('proteinAlignment').innerText =  protSimTable.row(this).data()[3].replace(/ /g, '\u00a0') 
        }
	} );
    
    $('#testLigandDataTable tbody').on('click', 'tr', function () {
        if (testLigandDataTable.row(this).length > 0) { 
            var molRegNo = testLigandDataTable.row(this).data()[6];
            document.getElementById("testLigandMolRegNo").innerText = molRegNo;
            $.ajax({
                url: "http://hadoop1:5432/getSmilesSVG/" + molRegNo + "/mol",
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
            var compId = testProteinDataTable.row(this).data()[6];
            var accession = testProteinDataTable.row(this).data()[7];
            document.getElementById("selectedTestProtein").innerText = accession;

            $.ajax({
            url: "http://hadoop1:5432/getProteinTestBindings/" + compId,
            type: "get",
            datatype: "json",
            success: function(response) {  
                var testProteinBindings = document.getElementById("testProteinBindings")              
                testProteinBindings.innerHTML = "<ul>"
                targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    testProteinBindings.innerHTML += ("<li><samp>" + targets[i] + "</samp></li>");
                }   
                testProteinBindings.innerHTML += ("</ul>");
                document.getElementById("testProteinBindingsTotal").innerText = "( " + targets.length + " )";
            }
            });
        }
	} );
    
    $('#ligandExperimentButton').click(function () {
        var molRegNo = document.getElementById("testLigandMolRegNo").innerText;
        document.getElementById("queryLigandMolRegNo").innerText = molRegNo;
        var fp = document.getElementById("ligFP1").value;
        var sim = document.getElementById("ligSim1").value;
        var th = document.getElementById("ligTH1").value;
        $.ajax({
            url: "http://hadoop1:5432/getSmilesSVG/" + molRegNo + "/mol",
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
        url: "http://hadoop1:5432/doLigandExperiment/" + molRegNo + "/" + fp + "/" + sim + "/" + th,
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
            $('#tabs').tabs({ active: 3 });
        }
        });        
	} );

    $('#ligSmilesExpRender').click(function () {
        var smiles = document.getElementById("ligSmilesExp").value;
        $.ajax({
            url: "http://hadoop1:5432/getSmilesSVG/" + smiles + "/smiles",
            type: "get",
            datatype: "json",
            success: function(response) {
                document.getElementById('ligSmilesSvg').innerHTML = response                                    
            }
        });      
	} );

    $('#ligSmilesExpRun').click(function () {
        var smiles = document.getElementById("ligSmilesExp").value;
        var fp = document.getElementById("ligFP1").value;
        var sim = document.getElementById("ligSim1").value;
        var th = document.getElementById("ligTH1").value;

        // render SVG
        $.ajax({
            url: "http://hadoop1:5432/getSmilesSVG/" + smiles + "/smiles",
            type: "get",
            datatype: "json",
            success: function(response) {
                document.getElementById('ligSmilesSvg').innerHTML = response;
                document.getElementById('queryLigandSVG').innerHTML = response;
            }
        });

        // check if ligand is in chembl
        $.ajax({
            url: "http://hadoop1:5432/isLigandInChEMBL/" + smiles,
            type: "get",
            datatype: "json",
            success: function(response) {
                document.getElementById("queryLigandMolRegNo").innerText = "In ChEMBL: " + response;                
            }
        });    

        // run experiment using smiles
        $.ajax({
        url: "http://hadoop1:5432/doLigandExperimentFromSmiles/" + smiles + "/" + fp + "/" + sim + "/" + th,
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
            $('#tabs').tabs({ active: 3 });
        }
        });
	} );

    $('#testProteinDataTable tbody').on('dblclick', 'tr', function () {
        if (testProteinDataTable.row(this).length > 0) { 
            var compId = testProteinDataTable.row(this).data()[6];
            var accession = testProteinDataTable.row(this).data()[7];
            document.getElementById("queryTestProtein").innerText = accession;

            $.ajax({
            url: "http://hadoop1:5432/getProteinTestBindings/" + compId,
            type: "get",
            datatype: "json",
            success: function(response) {  
                var testProteinBindings = document.getElementById("queryProteinBindings")              
                testProteinBindings.innerHTML = "<ul>"
                targets = JSON.parse(response)
                for (var i = 0; i < targets.length; i++) {
                    testProteinBindings.innerHTML += ("<li><samp>" + targets[i] + "</samp></li>");
                }   
                testProteinBindings.innerHTML += ("</ul>");
                document.getElementById("queryProteinBindingsTotal").innerText = "( " + targets.length + " )";
            }
            });

            $.ajax({
			url: "http://hadoop1:5432/doProteinExperiment/" + compId,
			type: "get",
			datatype: "json",			
			success: function (response) {				
				var res = JSON.parse(response);
                protSimTable.rows().clear();
                protSimTable.column(3).visible(false);
                protSimTable.rows.add(res).draw();

                var knownProteinsUniqueMolRegNo = document.getElementById("queryKnownProteinBindings");
                knownProteinsUniqueMolRegNo.innerHTML = "<ul>";
                var uniqueMolRegNo = protSimTable.column(5).data().unique().toArray().sort(sortNumber);
                for (var i = 0; i < uniqueMolRegNo.length; i++) {
                    knownProteinsUniqueMolRegNo.innerHTML += ("<li><samp>" + uniqueMolRegNo[i] + "</samp></li>");
                }
                knownProteinsUniqueMolRegNo.innerHTML += ("</ul>");
                document.getElementById("queryKnownProteinBindingsTotal").innerText = "( " + uniqueMolRegNo.length + " )";

                $('#tabs').tabs({ active: 5 });
			}
		    });
        }
	} );
    
    // http://gasparesganga.com/labs/jquery-loading-overlay/
    /*$(document).ajaxStart(function(){
        $.LoadingOverlay("show");
    });

    $(document).ajaxStop(function(){
        $.LoadingOverlay("hide");
    });*/
});

