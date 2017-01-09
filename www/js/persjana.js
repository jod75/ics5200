    var _map;

    function initialize() {
        var mapProp = {
            center: new google.maps.LatLng(35.9, 14.4),
            zoom: 12,
            mapTypeId: google.maps.MapTypeId.ROADMAP
        };
        
        var map = new google.maps.Map(document.getElementById("googleMap"), mapProp);           

        google.maps.event.addListener(map, 'mousemove', function (event) {
            var lat = event.latLng.lat();            
            var lng = event.latLng.lng();            
            $("#lat").text(lat.toFixed(4));
            $("#lon").text(lng.toFixed(4));        
        });

        google.maps.event.addListener(map, 'click', function (event) {
            showSearchText();
            showLatLngDetails(event.latLng.lat(), event.latLng.lng());
        })

        _map = map;
    }

    jQuery(function ($) {
        $('#panels').enhsplitter({ position: '85%' });
    });

    $(document).ready(function () {
        // on document ready, get the callsigns recorded (last 24hrs)
        $.ajax({
            url: "./api/persjana.py",
            type: "post",
            datatype: "json",
            data: { 'key1': 'getCallsigns'},
            success: function (response) {
                callsigns = JSON.parse(response.data)
                // fill drop down
                var dropdown = $("#callsigns");
                dropdown.empty(); 
                for (var i = 0; i < callsigns.length; i++) {
                    dropdown.append("<option>" + callsigns[i]['_id']['from'] + "</option>");
                }          
            }
        });

        $("#locateCallsign").click(function () {
            showSearchText();
            $.ajax({
                url: "./api/persjana.py",
                type: "post",
                datatype: "json",
                data: { 'key1': 'getCallsignsLocation', 'callsign': $("#callsigns").val() },
                success: function (response) {
                    detailsoflocation = JSON.parse(response.data)
                    if (detailsoflocation.length > 0) {
                        lat = detailsoflocation[0]['latitude'];
                        lng = detailsoflocation[0]['longitude'];
                        _map.setCenter(new google.maps.LatLng(lat, lng));
                        showLatLngDetails(lat, lng);
                    }
                }
            });
        });

        $("#allActivity").click(function () {
            $.ajax({
                url: "./api/persjana.py",
                type: "post",
                datatype: "json",
                data: {
                    'key1': 'getActivity',
                    'longitude1': _map.getBounds().getNorthEast().lng(),
                    'latitude1': _map.getBounds().getNorthEast().lat(),
                    'longitude2': _map.getBounds().getSouthWest().lng(),
                    'latitude2': _map.getBounds().getSouthWest().lat()
                },
                success: function (response) {
                    showMapData(response.data, '#FF0000')
                }
            });
        });

        $("#aprsRepeaters").click(function () {
            $.ajax({
                url: "./api/persjana.py",
                type: "post",
                datatype: "json",
                data: {
                    'key1': 'getActivityBySymbol',
                    'longitude1': _map.getBounds().getNorthEast().lng(),
                    'latitude1': _map.getBounds().getNorthEast().lat(),
                    'longitude2': _map.getBounds().getSouthWest().lng(),
                    'latitude2': _map.getBounds().getSouthWest().lat(),
                    'symbol': 'r',
                    'symbol_table': '/'
                },
                success: function (response) {
                    showMapData(response.data, '#000000')
                }
            });
        });

        $("#citizenWX").click(function () {
            $.ajax({
                url: "./api/persjana.py",
                type: "post",
                datatype: "json",
                data: {
                    'key1': 'getActivityBySymbol',
                    'longitude1': _map.getBounds().getNorthEast().lng(),
                    'latitude1': _map.getBounds().getNorthEast().lat(),
                    'longitude2': _map.getBounds().getSouthWest().lng(),
                    'latitude2': _map.getBounds().getSouthWest().lat(),
                    'symbol': '_',
                    'symbol_table': '/'
                },
                success: function (response) {
                    showMapData(response.data, '#0000FF')
                }
            });
        });

	$("#noaaStations").click(function () {
            $.ajax({
                url: "./api/persjana.py",
                type: "post",
                datatype: "json",
                data: {
                    'key1': 'getNOAAStations',
                    'longitude1': _map.getBounds().getNorthEast().lng(),
                    'latitude1': _map.getBounds().getNorthEast().lat(),
                    'longitude2': _map.getBounds().getSouthWest().lng(),
                    'latitude2': _map.getBounds().getSouthWest().lat()
                },
                success: function (response) {
                    showMapNOAAData(response.data, '#FF00FF')
                }
            });
        });


    });

    // Google maps events
    google.maps.event.addDomListener(window, 'load', initialize);

    function showMapData(data, c) {
        for (var act in data) {
            // Add the circle for this city to the map.
            var cityCircle = new google.maps.Circle({
                strokeColor: c,
                strokeOpacity: 0.8,
                strokeWeight: 2,
                fillColor: c,
                fillOpacity: 0.35,
                map: _map,
                center: { lat: data[act]._id.lat, lng: data[act]._id.lon },
                radius: Math.sqrt(data[act].activity) * 100
            });
        }
    }

    function showMapNOAAData(data, c) {
        for (var act in data) {
            // Add the circle for this city to the map.
            var cityCircle = new google.maps.Circle({
                strokeColor: c,
                strokeOpacity: 0.8,
                strokeWeight: 2,
                fillColor: c,
                fillOpacity: 0.35,
                map: _map,
                center: { lat: data[act]._id.latitude, lng: data[act]._id.longitude },
                radius: Math.sqrt(100) * 100
            });
        }
    }

    function showLatLngDetails(lat, lng) {
        $.ajax({
            url: "./api/persjana.py",
            type: "post",
            datatype: "json",
            data: { 'key1': 'getPoi', 'latitude': lat, 'longitude': lng, 'radius': $("#radius").val() },
            success: function (response) {
                // display POI
                poidata = unescape(JSON.parse(response.data))
                $("#square_text").text(poidata)
            },
            error: function () {
                $("#square_text").text('Failed to get data')
            }
        })
        $.ajax({
            url: "./api/persjana.py",
            type: "post",
            datatype: "json",
            data: { 'key1': 'getCurrentWeather', 'latitude': lat, 'longitude': lng },
            success: function (response) {
                // display Weather
                weatherdata = JSON.parse(response.data)
                $("#circle_text").text(weatherdata)
            },
            error: function () {
                $("#circle_text").text('Failed to get data')
            }
        })
        $.ajax({
            url: "./api/persjana.py",
            type: "post",
            datatype: "json",
            data: { 'key1': 'getCurrentLocationInfo', 'latitude': lat, 'longitude': lng, 'radius': $("#radius").val()},
            success: function (response) {
                // display Location Info
                locationinfodata = unescape(JSON.parse(response.data))
                $("#triangle_text").text(locationinfodata )
            },
            error: function () {
                $("#triangle_text").text('Failed to get data')
            }
        })
    }

    function showSearchText() {
        $('#circle_text').text('Fetching...')
        $('#triangle_text').text('Searching for data...')
        $('#square_text').text('Please wait...')
    }