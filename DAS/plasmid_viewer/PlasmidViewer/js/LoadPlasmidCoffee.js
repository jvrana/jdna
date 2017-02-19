// Generated by CoffeeScript 1.12.2
(function() {
  'use strict';
  var app, assign_features, feature_loc, get_feature_length, sanitize_features, sort_features;

  app = angular.module('myApp', ['angularplasmid']);

  feature_loc = function(feature, start, end, length) {
    var e, s;
    s = start(feature);
    e = end(feature);
    console.log(s, length, e);
    if (s > length) {
      s = s - length;
    }
    if (e > length) {
      e = e - length;
    }
    feature.plasmid_start = s;
    feature.plasmid_end = e;
    return [s, e];
  };

  sanitize_features = function(features, start, end, length) {
    var f, i, len, results;
    results = [];
    for (i = 0, len = features.length; i < len; i++) {
      f = features[i];
      results.push(feature_loc(f, start, end, length));
    }
    return results;
  };

  get_feature_length = function(a, length) {
    var x;
    x = a.plasmid_end - a.plasmid_start;
    if (x < 0) {
      x = x + length;
    }
    if (x === (0/0)) {
      console.log('NNNN', a.plasmid_start, a.plasmid_end);
    }
    return x;
  };

  sort_features = function(features, length) {
    return features.sort(function(a, b) {
      return get_feature_length(a, length) - get_feature_length(b, length);
    }).reverse();
  };

  assign_features = function(features, length) {
    var af, assigned_features, feature, feature_overlaps, i, j, len, len1, max_feature_shell;
    feature_overlaps = function(f1, f2, length) {
      var e, origin, overlap, ref, ref1, ref2, ref3, ref4, ref5, s;
      origin = false;
      overlap = false;
      s = f2.plasmid_start;
      e = f2.plasmid_end;
      if (f1.shell === f2.shell) {
        if (s > e) {
          origin = true;
        }
        if (!origin) {
          if ((s <= (ref = f1.plasmid_start) && ref < e)) {
            overlap = true;
          }
          if ((s < (ref1 = f1.plasmid_end) && ref1 <= e)) {
            overlap = true;
          }
        } else {
          if ((s <= (ref2 = f1.plasmid_start) && ref2 <= length) || (0 < (ref3 = f1.plasmid_start) && ref3 < e)) {
            overlap = true;
          }
          if ((s < (ref4 = f1.plasmid_end) && ref4 <= length) || (0 <= (ref5 = f1.plasmid_end) && ref5 <= e)) {
            overlap = true;
          }
        }
        if (overlap) {
          console.log('\t' + f1.name + ' (' + f1.plasmid_start + '-' + f1.plasmid_end + ') overlaps with ' + f2.name + ' (' + f2.plasmid_start + '-' + f2.plasmid_end + ') in shell ' + f1.shell);
        }
      }
      return overlap;
    };
    max_feature_shell = 0;
    assigned_features = [];
    for (i = 0, len = features.length; i < len; i++) {
      feature = features[i];
      if (!feature.shell) {
        feature.shell = 1;
        for (j = 0, len1 = assigned_features.length; j < len1; j++) {
          af = assigned_features[j];
          if (feature_overlaps(feature, af) || feature_overlaps(af, feature)) {
            feature.shell = feature.shell + 1;
          }
        }
        if (feature.shell > max_feature_shell) {
          max_feature_shell = feature.shell;
        }
        console.log("Assigning " + feature.name + ' (' + feature.plasmid_start + '-' + feature.plasmid_end + ' to ' + feature.shell);
        assigned_features.push(feature);
      }
    }
    return max_feature_shell;
  };

  app.controller('LoadPlasmids', function($scope, $http) {
    $http({
      method: 'POST',
      url: 'js/plasmiddata.json'
    }).success(function(data) {
      var end, f, i, len, padding, ref, start;
      $scope.plasmid = data;
      $scope.plasmid_size = 150;
      $scope.feature_width = 10;
      $scope.feature_padding = 3;
      start = function(x) {
        return x.start;
      };
      end = function(x) {
        return x.stop;
      };
      sanitize_features($scope.plasmid.features, start, end, data.length);
      console.log("Sorting features");
      data.features = sort_features(data.features, data.length);
      ref = data.features;
      for (i = 0, len = ref.length; i < len; i++) {
        f = ref[i];
        console.log("L " + get_feature_length(f, data.length), f.plasmid_start, f.plasmid_end, data.length);
      }
      $scope.max_feature_shell = assign_features($scope.plasmid.features, data.length);
      padding = 30;
      return $scope.contig_track_size = padding + $scope.plasmid_size + $scope.max_feature_shell * ($scope.feature_width + $scope.feature_padding);
    });
    $http({
      method: 'POST',
      url: '../../alignment_viewer/data.json'
    }).success(function(data) {
      var end, start;
      $scope.contigs = data;
      $scope.contig_padding = 10;
      start = function(x) {
        return x.query.start;
      };
      end = function(x) {
        return x.query.end;
      };
      sanitize_features(data.contigs, start, end, $scope.plasmid.length);
      data.contigs = sort_features(data.features, data.length);
      return $scope.max_contig_shell = assign_features(data.contigs, data.length);
    });
    $http({
      method: 'POST',
      url: '../../alignment_viewer/primer_data.json'
    }).success(function(data) {
      return $scope.primers = data;
    });
    $http({
      method: 'POST',
      url: '../../alignment_viewer/bestassembly.json'
    }).success(function(data) {
      return $scope.assembly = data;
    });
    return $scope.color = function() {
      var c;
      c = "rgba(" + (Math.ceil(Math.random() * 255)) + "," + (Math.ceil(Math.random() * 255)) + "," + (Math.ceil(Math.random() * 255)) + "," + 0.9 + ")";
      return c;
    };
  });

}).call(this);

//# sourceMappingURL=LoadPlasmidCoffee.js.map
