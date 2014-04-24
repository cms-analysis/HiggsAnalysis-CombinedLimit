(function () {
    'use strict';

    var url = '/datacards';

    angular.module('datacardsApp', [
        'blueimp.fileupload'
    ])
        .config([
            '$httpProvider', 'fileUploadProvider',
            function ($httpProvider, fileUploadProvider) {
                delete $httpProvider.defaults.headers.common['X-Requested-With'];
                fileUploadProvider.defaults.redirect = window.location.href.replace(
                    /\/[^\/]*$/,
                    '/cors/result.html?%s'
                );
                angular.extend(fileUploadProvider.defaults, {
                    maxFileSize: 50000000,
                    acceptFileTypes: /(\.|\/)(txt|root)$/i
                });
            }
        ])

        .controller('DemoFileUploadController', [
            '$scope', '$http', '$filter', '$window',
            function ($scope, $http) {
                $scope.options = {
                    url: url
                };
                $scope.loadingFiles = true;
                $http.get(url)
                    .then(
                        function (response) {
                            $scope.loadingFiles = false;
                            $scope.queue = response.data.files || [];
                        },
                        function () {
                            $scope.loadingFiles = false;
                        }
                    );
            }
        ])

        .controller('FileDestroyController', [
            '$scope', '$http',
            function ($scope, $http) {
                var file = $scope.file,
                    state;
                if (file.url) {
                    file.$state = function () {
                        return state;
                    };
                    file.$destroy = function () {
                        state = 'pending';
                        return $http({
                            url: file.deleteUrl,
                            method: file.deleteType
                        }).then(
                            function () {
                                state = 'resolved';
                                $scope.clear(file);
                            },
                            function () {
                                state = 'rejected';
                            }
                        );
                    };
                } else if (!file.$cancel && !file._index) {
                    file.$cancel = function () {
                        $scope.clear(file);
                    };
                }
            }
        ])

        .controller('ViewDatacardsController', [
            '$scope', '$http', '$filter', '$window', '$parse',
            function ($scope, $http) {
                $scope.refresh = function(){
                    $scope.datacards = [];
                    $scope.total = 0;
                    $http.get(url+"/list")
                    .success(function(data) {
                    if ( angular.isArray(data.files) ) {
                        $scope.datacards = data.files;
                    }
                    else {
                        $scope.datacards = [data.files];
                    }
                    $scope.total = $scope.datacards.length;
                    })
                };
                $scope.get_datacard = function(filename){
                    var lines = [];
                    $http.get(url+'/'+filename)
                    .success(function(data) {
                        show_datacard(data);
                    })
                };
            }
        ]);
}());