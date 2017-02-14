/**
 * Created by Justin on 2/13/17.
 */

'use strict';

// Declare app level module which depends on views, and components
var app = angular.module('myApp', ['angularplasmid']);

app.controller('LoadPlasmids', function($scope, $http)
{
    $http({method: 'POST', url: 'js/data.json'}).success(function(data)
    {
        $scope.plasmids = data; // response data
    });
});