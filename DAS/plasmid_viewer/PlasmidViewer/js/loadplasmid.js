/**
 * Created by Justin on 2/13/17.
 */

'use strict';

// Declare app level module which depends on views, and components
var app = angular.module('myApp', ['angularplasmid']);

app.controller('LoadPlasmids', function($scope, $http)
{
    $http({method: 'POST', url: 'js/plasmiddata.json'}).success(function(data)
    {
        $scope.plasmids = data; // response data
    });
    $http({method: 'POST', url: '../../alignment_viewer/data.json'}).success(function(data)
    {
        $scope.contigs = data; // response data
    });

    $http({method: 'POST', url: '../../alignment_viewer/primer_data.json'}).success(function(data)
    {
        $scope.primers = data; // response data
    });

    $http({method: 'POST', url: '../../alignment_viewer/bestassembly.json'}).success(function(data)
    {
        $scope.assembly = data; // response data
    });

    $scope.color = function() {
        return "rgba(" + (Math.ceil(Math.random() * 255)) + "," +
            (Math.ceil(Math.random() * 255)) + "," +
            (Math.ceil(Math.random() * 255)) + "," +
            0.5 + ")"
    }

    $scope.getNumber = function() {
        $scope.num = (Math.ceil(Math.random() * 255));
    };
});
